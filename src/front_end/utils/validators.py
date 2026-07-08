import re

def validate_ncbi_accession(accession: str) -> bool:
    """
    Validates if a given string matches the standard NCBI Accession format.
    E.g. CP031675.1, NM_001256799.3, ABCDE123456
    """
    accession = accession.strip()
    if not accession:
        return False
        
    # Standard RefSeq and GenBank format (1-4 letters, optional underscore, followed by 5+ digits, optional version .1)
    # Examples:
    # U12345
    # AF123456
    # NC_000913.3
    # NZ_CP031675.1
    
    pattern = r"^[A-Z]{1,4}_?\d{5,}(\.\d+)?$"
    
    # Check for prefix (e.g. NZ_ or NC_) before the main pattern
    if accession.startswith("NZ_") or accession.startswith("NC_") or accession.startswith("NM_") or accession.startswith("NP_"):
        main_part = accession[3:]
        return bool(re.match(r"^[A-Z]{0,4}\d{5,}(\.\d+)?$", main_part))
        
    return bool(re.match(pattern, accession))

def parse_comma_separated_accessions(accessions_str: str) -> list:
    """
    Parses a comma-separated string of accessions and returns a list of valid ones.
    Returns a tuple of (valid_list, invalid_list).
    """
    parts = [p.strip() for p in accessions_str.split(',') if p.strip()]
    valid = []
    invalid = []
    for p in parts:
        if validate_ncbi_accession(p):
            valid.append(p)
        else:
            invalid.append(p)
    return valid, invalid
