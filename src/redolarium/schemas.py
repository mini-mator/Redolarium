# type: ignore
# Pydantic v2 schemas and validation models for GenoMeta v2.0
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Tuple, Literal
import json

SCHEMA_VERSION = "2.0.0"

class GenomeMetadata(BaseModel):
    genome_id: str = Field(..., description="Unique genome identifier")
    raw_fasta_path: str = Field(..., description="Path to the raw FASTA/GenBank assembly file")
    contig_count: int = Field(..., gt=0, description="Total number of contigs in the assembly")
    total_length_bp: int = Field(..., gt=0, description="Total nucleotide length of the genome in bp")
    execution_seed: int = Field(42, ge=0, description="Random seed used for deterministic execution tracking")

class SyntenyBlock(BaseModel):
    query_contig: str = Field(..., description="Query contig identifier")
    target_accession: str = Field(..., description="Accession of target reference strain")
    identity_percentage: float = Field(..., ge=0.0, le=100.0, description="Percent identity across the syntenic region")
    jaccard_index: float = Field(..., ge=0.0, le=1.0, description="Jaccard similarity index of orthology clusters")

class BGCOutput(BaseModel):
    genome_id: str = Field(..., description="Source genome ID")
    bgc_id: str = Field(..., description="Biosynthetic Gene Cluster ID")
    predicted_type: str = Field(..., description="Predicted type of BGC (e.g., NRPS, PKS, Lantipeptide)")
    start_coord: int = Field(..., ge=0, description="BGC starting nucleotide coordinate")
    end_coord: int = Field(..., ge=0, description="BGC ending nucleotide coordinate")
    pfam_domains: List[str] = Field(default_factory=list, description="List of detected Pfam domains in the cluster")
    mibig_hits: List[Dict[str, Any]] = Field(default_factory=list, description="List of MIBiG similarity match results")

class HGTOutput(BaseModel):
    genome_id: str = Field(..., description="Source genome ID")
    contig_id: str = Field(..., description="Contig identifier")
    coordinates: Tuple[int, int] = Field(..., description="Sliding window coordinates tuple (start, end)")
    gc_deviation: float = Field(..., description="GC content deviation from global genome average")
    tnf_distance: float = Field(..., ge=0.0, description="Euclidean distance of tetranucleotide frequencies")
    codon_bias_index: float = Field(..., ge=0.0, le=1.0, description="Calculated codon usage bias index")
    mge_proximity_bp: int = Field(..., ge=0, description="Distance to the closest MGE / prophage element in bp")
    confidence_assignment: Literal['High', 'Medium', 'Low'] = Field(..., description="Composite HGT confidence level")

class DockingOutput(BaseModel):
    genome_id: str = Field(..., description="Source genome ID")
    target_protein_id: str = Field(..., description="Receptor peptidase protein ID")
    ligand_id: str = Field(..., description="Ligand/precursor peptide BGC identifier")
    binding_affinity_kcal_mol: float = Field(..., description="Top mode binding affinity in kcal/mol (negative is stronger)")
    contact_residues: str = Field(..., description="List of key interaction contact residues")
    engine_used: str = Field(..., description="Engine used for prediction (e.g., smina, vinardo, none_alphafold_fallback)")
    execution_status: str = Field(..., description="Execution status (e.g., Success, Degraded_Coordinates_Only)")


# ==========================================
# STRICT UPSTREAM PARSING & SERIALIZATION EXAMPLE
# ==========================================
def parse_and_serialize_hgt_raw(raw_data: Dict[str, Any]) -> str:
    """
    Parses a raw output dictionary (from an upstream HGT tool parser) 
    into a type-validated Pydantic model and serializes it to JSON,
    intercepting any malformed inputs.
    """
    # Upstream data could contain raw values; we enforce type conversion and validate bounds
    validated_model = HGTOutput(
        genome_id=str(raw_data.get("genome_id")),
        contig_id=str(raw_data.get("contig_id")),
        coordinates=(int(raw_data.get("start", 0)), int(raw_data.get("end", 0))),
        gc_deviation=float(raw_data.get("gc_deviation", 0.0)),
        tnf_distance=float(raw_data.get("tnf_distance", 0.0)),
        codon_bias_index=float(raw_data.get("codon_bias_index", 0.0)),
        mge_proximity_bp=int(raw_data.get("mge_proximity_bp", 999999)),
        confidence_assignment=raw_data.get("confidence_assignment", "Low")
    )
    return validated_model.model_dump_json(indent=2)
