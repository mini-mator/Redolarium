import os
import re
import csv
import time
from Bio.Blast import NCBIXML
import io
from Bio import Entrez

# Reference: Konstantinidis & Tiedje 2005 (doi:10.1073/pnas.0409727102)
# An 85% amino acid identity threshold on core housekeeping marker genes ensures 
# genus/species-level boundary mapping, preventing taxonomic spillover (Rule R4).
IDENTITY_THRESHOLD = 85.0

def sweep_local_db(program, query_seq, query_record, logger):
    """
    Tier 1 Sweep: Run BLAST against the local ref_markers_db.
    Parse alignments using regex to extract local reference accessions.
    Enforce strict % Identity Threshold Gate.
    """
    from redolarium.annotation import run_blast_query
    
    local_refs = []
    seen_ignored = set()
    blast_record = None
    
    # BUGFIX: ref_markers_db is a protein database. If program is blastn, this will crash.
    if program == "blastn":
        logger.info("[Local_DB] Bypassing Tier 1 sweep. Cannot run blastn against protein database (ref_markers_db).")
        return local_refs, blast_record
        
    try:
        xml_str = run_blast_query(program, "ref_markers_db", query_seq, logger)
        if xml_str:
            blast_record = NCBIXML.read(io.StringIO(xml_str))
            for alignment in blast_record.alignments:
                title = alignment.title
                identity_pct = 0.0
                if alignment.hsps:
                    hsp = alignment.hsps[0]
                    identities = hsp.identities
                    align_len = hsp.align_length
                    identity_pct = (identities / align_len) * 100.0 if align_len > 0 else 0.0
                    
                match_local = re.search(r'([a-zA-Z0-9\.]+)_([a-zA-Z0-9]+)\s+\[([^\]]+)\]', title, re.IGNORECASE)
                if match_local:
                    acc_val = match_local.group(1)
                    species_val = match_local.group(3)
                    ref_name = f"{species_val} ({acc_val})"
                    if identity_pct >= IDENTITY_THRESHOLD:
                        if ref_name not in local_refs and acc_val != query_record.id:
                            local_refs.append(ref_name)
                            logger.info(f"[Local_DB] Valid Sweep hit found: {species_val} (ID: {identity_pct:.1f}%).")
                    else:
                        if ref_name not in seen_ignored:
                            seen_ignored.add(ref_name)
                            logger.info(f"[Local_DB] Ignored low-identity hit: {species_val} (ID: {identity_pct:.1f}%).")
    except Exception as e:
        logger.warning(f"Local sweep failed: {e}")
        
    return local_refs, blast_record

class FileLock:
    """A cross-process file-based lock to serialize remote NCBI queries."""
    def __init__(self, lock_path):
        self.lock_path = lock_path

    def acquire(self):
        while True:
            try:
                fd = os.open(self.lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.close(fd)
                return
            except FileExistsError:
                time.sleep(1.5)

    def release(self):
        try:
            os.remove(self.lock_path)
        except OSError:
            pass

def get_nucleotide_acc_from_assembly(assembly_id):
    """Link NCBI Assembly ID to Nuccore ID and return AccessionVersion."""
    try:
        link_handle = Entrez.elink(dbfrom="assembly", db="nuccore", id=assembly_id)
        link_record = Entrez.read(link_handle)
        links = []
        try:
            link_set = link_record[0]["LinkSetDb"]
            for db_link in link_set:
                if db_link["DbTo"] == "nuccore":
                    links = [link["Id"] for link in db_link["Link"]]
                    break
        except Exception:
            pass
        if links:
            sum_handle = Entrez.esummary(db="nuccore", id=links[0])
            record_sum = Entrez.read(sum_handle)
            acc = record_sum[0].get('AccessionVersion', record_sum[0].get('Caption', ''))
            if acc:
                return acc
    except Exception:
        pass
    return None

def run_remote_blast_expansion(program, query_seq, current_references, limit, query_record, logger, email="researcher@example.com"):
    """
    Tier 2 Sweep: If local sweep yielded < 50 references, query remote BLAST.
    Parse remote hits incrementally to find unique references up to 50.
    Enforce strict % Identity Threshold Gate and exhaustion flags.
    """
    from Bio.Blast import NCBIWWW
    import time
    
    remote_refs = []
    search_space_exhausted = False
    
    # Establish lock path in the resources directory
    workspace_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    lock_dir = os.path.join(workspace_dir, "resources")
    os.makedirs(lock_dir, exist_ok=True)
    lock_path = os.path.join(lock_dir, "ncbi_blast.lock")
    lock = FileLock(lock_path)
    
    try:
        logger.info("Waiting for lock to perform remote NCBI BLAST query (preventing parallel IP throttling)...")
        lock.acquire()
        logger.info("Lock acquired. Executing remote NCBI BLAST query for reference expansion...")
        result_handle = NCBIWWW.qblast(program, "nr" if program == "blastp" else "nt", query_seq, hitlist_size=100)
        xml_str = result_handle.read()
    finally:
        lock.release()
        
    try:
        
        if xml_str:
            blast_record = NCBIXML.read(io.StringIO(xml_str))
            
            if len(blast_record.alignments) == 0:
                search_space_exhausted = True
                logger.warning("Remote BLAST returned 0 alignments. Search space exhausted.")
                return remote_refs, search_space_exhausted
                
            idx = 0
            total_alignments = len(blast_record.alignments)
            
            while len(current_references) + len(remote_refs) < limit:
                if idx >= total_alignments:
                    search_space_exhausted = True
                    break
                    
                alignment = blast_record.alignments[idx]
                accession = alignment.accession
                title = alignment.title
                idx += 1
                
                identity_pct = 0.0
                if alignment.hsps:
                    hsp = alignment.hsps[0]
                    identities = hsp.identities
                    align_len = hsp.align_length
                    identity_pct = (identities / align_len) * 100.0 if align_len > 0 else 0.0
                    
                match = re.search(r'\[([^\]]+)\]', title)
                if match:
                    species = match.group(1)
                    if identity_pct >= IDENTITY_THRESHOLD:
                        logger.info(f"[Remote_DB] Valid Sweep hit found: {species} (ID: {identity_pct:.1f}%).")
                        try:
                            Entrez.email = email
                            term = f"{species}[Organism] AND (complete genome[Assembly Level] OR chromosome[Assembly Level])"
                            handle = Entrez.esearch(db="assembly", term=term, retmax=limit)
                            record = Entrez.read(handle)
                            id_list = record["IdList"]
                            if id_list:
                                sum_handle = Entrez.esummary(db="assembly", id=",".join(id_list))
                                record_sum = Entrez.read(sum_handle)
                                doc_sums = record_sum['DocumentSummarySet']['DocumentSummary']
                                for doc_idx, doc in enumerate(doc_sums):
                                    level = doc.get('AssemblyStatus', doc.get('AssemblyLevel', ''))
                                    if level in ["Complete Genome", "Chromosome", "Scaffold", "Contig"]:
                                        name = doc.get('Organism', '')
                                        # Resolve true nucleotide accession
                                        assembly_id = id_list[doc_idx] if doc_idx < len(id_list) else None
                                        acc = get_nucleotide_acc_from_assembly(assembly_id) if assembly_id else None
                                        if not acc:
                                            acc = doc.get('AssemblyAccession', '')
                                        if name and acc:
                                            ref_name = f"{name} ({acc})"
                                            if (ref_name not in current_references and 
                                                ref_name not in remote_refs and 
                                                acc != query_record.id):
                                                remote_refs.append(ref_name)
                                                logger.info(f"[Remote_Fetched] Resolved reference: {ref_name}")
                                                break
                        except Exception as ex:
                            logger.warning(f"Failed to query assembly for species {species}: {ex}")
                    else:
                        logger.info(f"[Remote_DB] Ignored low-identity hit: {species} (ID: {identity_pct:.1f}%).")
                        
            if idx >= total_alignments:
                search_space_exhausted = True
                
    except Exception as e:
        logger.warning(f"Remote BLAST expansion failed: {e}")
        search_space_exhausted = True
        
    return remote_refs, search_space_exhausted

def export_provenance_report(out_dir, local_refs, remote_refs):
    """
    Write results_provenance.csv to out_dir, categorizing references into
    [Local_DB, Remote_Fetched] columns.
    """
    csv_path = os.path.join(out_dir, "results_provenance.csv")
    max_len = max(len(local_refs), len(remote_refs))
    
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Local_DB", "Remote_Fetched"])
        for i in range(max_len):
            loc = local_refs[i] if i < len(local_refs) else ""
            rem = remote_refs[i] if i < len(remote_refs) else ""
            writer.writerow([loc, rem])
