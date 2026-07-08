import os
import sys
import platform
import time

class PredictionResult:
    """
    Standardized result object for all Redolarium analysis modules.
    Tracks predictions, continuous confidence scores (0.0 to 1.0),
    provenance details, databases, execution parameters, and citations.
    """
    def __init__(self, prediction, confidence_score, algorithm, algorithm_version,
                 genome_closure_status="Closed", fallback_active=False, fallback_multiplier=0.5,
                 uncertainty=None, database=None, database_version=None, download_date=None,
                 checksum=None, source_url=None, evidence=None, limitations=None,
                 citations=None, runtime=0.0, warnings=None, metadata=None):
        self.prediction = prediction
        self.genome_closure_status = str(genome_closure_status)
        self.fallback_active = bool(fallback_active)
        self.fallback_multiplier = float(fallback_multiplier)
        
        # Calculate base score
        raw_score = float(confidence_score)
        
        # Apply mandatory assembly fragmentation penalty if the genome is a draft
        if self.genome_closure_status.lower() == "draft":
            raw_score *= 0.80  # Mandatory 20% downrate for fragmented assemblies
            
        # Apply automatic fallback downrating if active
        if self.fallback_active:
            self.confidence_score = max(0.0, min(1.0, raw_score * self.fallback_multiplier))
        else:
            self.confidence_score = max(0.0, min(1.0, raw_score))
            
        self.confidence_label = self._derive_label(self.confidence_score)
        self.manual_intervention_required = (self.confidence_score < 0.3)
        self.uncertainty = uncertainty
        self.algorithm = algorithm
        self.algorithm_version = algorithm_version
        
        # Database Provenance (FAIR Principles)
        self.database = database
        self.database_version = database_version
        self.download_date = download_date
        self.checksum = checksum
        self.source_url = source_url
        
        self.evidence = evidence or []
        self.limitations = limitations or []  # These will map to "Systemic Uncertainties" in reporting
        if self.genome_closure_status.lower() == "draft":
            self.limitations.append("Input assembly is a draft (multi-contig). Downstream synteny and operon mapping may have reduced resolution.")
        if self.fallback_active:
            self.limitations.append(f"Tool executed in fallback mode. Confidence scaled by {self.fallback_multiplier}x.")
            
        self.citations = citations or []
        self.runtime = runtime
        self.warnings = warnings or []
        self.metadata = metadata or {}
        
        # Determine psutil availability for RAM measurement
        ram_gb = "N/A"
        try:
            import psutil
            ram_gb = round(psutil.virtual_memory().total / (1024.0**3), 2)
        except ImportError:
            pass

        # FAIR execution provenance data
        self.provenance = {
            "command_line": " ".join(sys.argv),
            "parameters": self.metadata.get("parameters", {}),
            "container_hash": os.environ.get("DOCKER_CONTAINER_HASH", "N/A"),
            "cpu_model": platform.processor() or platform.machine(),
            "ram_total_gb": ram_gb,
            "os_platform": platform.system(),
            "python_version": platform.python_version()
        }

    def _derive_label(self, score):
        if score >= 0.9: return "High"
        elif score >= 0.7: return "Moderate"
        elif score >= 0.5: return "Low"
        else: return "Very Low"

    def to_dict(self):
        return {
            "prediction": self.prediction,
            "confidence_score": self.confidence_score,
            "confidence_label": self.confidence_label,
            "manual_intervention_required": self.manual_intervention_required,
            "genome_closure_status": self.genome_closure_status,
            "fallback_active": self.fallback_active,
            "fallback_multiplier": self.fallback_multiplier,
            "uncertainty": self.uncertainty,
            "algorithm": self.algorithm,
            "algorithm_version": self.algorithm_version,
            "database": self.database,
            "database_version": self.database_version,
            "download_date": self.download_date,
            "checksum": self.checksum,
            "source_url": self.source_url,
            "evidence": self.evidence,
            "limitations": self.limitations,
            "citations": self.citations,
            "runtime": self.runtime,
            "warnings": self.warnings,
            "metadata": self.metadata,
            "provenance": self.provenance
        }

    def __iter__(self):
        if isinstance(self.prediction, list):
            return iter(self.prediction)
        elif isinstance(self.prediction, dict):
            return iter(self.prediction.items())
        raise TypeError(f"Prediction of type {type(self.prediction)} is not iterable")
