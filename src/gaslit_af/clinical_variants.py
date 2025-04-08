"""
Clinical Variants Module for GASLIT-AF Variant Analysis.

This module handles the processing, validation, and integration of clinical variant data
according to the defined JSON schema.
"""

import os
import json
import logging
import jsonschema
from pathlib import Path
from typing import Dict, List, Any, Optional, Union

# Configure logging
log = logging.getLogger("gaslit-af")

class ClinicalVariantManager:
    """Manager for clinical variant data processing and validation."""
    
    def __init__(self, schema_path: Optional[str] = None):
        """
        Initialize the clinical variant manager.
        
        Args:
            schema_path: Path to the JSON schema file. If None, uses the default schema.
        """
        if schema_path is None:
            # Use the default schema included with the package
            schema_path = Path(__file__).parent / "schemas" / "clinical_variants.json"
        
        self.schema_path = Path(schema_path)
        self.schema = self._load_schema()
        self.conditions_data = {"conditions": []}
        
    def _load_schema(self) -> Dict:
        """
        Load the JSON schema from file.
        
        Returns:
            Dict: The loaded JSON schema
        """
        try:
            with open(self.schema_path, 'r') as f:
                schema = json.load(f)
            log.info(f"Loaded clinical variants schema from {self.schema_path}")
            return schema
        except Exception as e:
            log.error(f"Error loading schema from {self.schema_path}: {e}")
            raise
    
    def load_conditions(self, conditions_file: Union[str, Path]) -> Dict:
        """
        Load conditions data from a JSON file.
        
        Args:
            conditions_file: Path to the conditions JSON file
            
        Returns:
            Dict: The loaded conditions data
        """
        try:
            with open(conditions_file, 'r') as f:
                self.conditions_data = json.load(f)
            
            # Validate against schema
            self.validate_conditions(self.conditions_data)
            
            log.info(f"Loaded {len(self.conditions_data.get('conditions', []))} conditions from {conditions_file}")
            return self.conditions_data
        except Exception as e:
            log.error(f"Error loading conditions from {conditions_file}: {e}")
            raise
    
    def validate_conditions(self, conditions_data: Dict) -> bool:
        """
        Validate conditions data against the JSON schema.
        
        Args:
            conditions_data: The conditions data to validate
            
        Returns:
            bool: True if validation passes, raises exception otherwise
        """
        try:
            jsonschema.validate(instance=conditions_data, schema=self.schema)
            log.info("Conditions data validated successfully against schema")
            return True
        except jsonschema.exceptions.ValidationError as e:
            log.error(f"Validation error: {e}")
            raise
    
    def add_condition(self, condition: Dict) -> None:
        """
        Add a new condition to the conditions data.
        
        Args:
            condition: The condition data to add
        """
        # Validate the single condition against the Condition schema
        condition_schema = self.schema["definitions"]["Condition"]
        try:
            jsonschema.validate(instance=condition, schema=condition_schema)
            self.conditions_data["conditions"].append(condition)
            log.info(f"Added condition: {condition.get('name')}")
        except jsonschema.exceptions.ValidationError as e:
            log.error(f"Validation error for condition {condition.get('name', 'unknown')}: {e}")
            raise
    
    def save_conditions(self, output_file: Union[str, Path]) -> None:
        """
        Save the current conditions data to a JSON file.
        
        Args:
            output_file: Path to save the conditions data
        """
        try:
            # Ensure the directory exists
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_path, 'w') as f:
                json.dump(self.conditions_data, f, indent=2)
            
            log.info(f"Saved {len(self.conditions_data.get('conditions', []))} conditions to {output_file}")
        except Exception as e:
            log.error(f"Error saving conditions to {output_file}: {e}")
            raise
    
    def get_condition_by_gene(self, gene: str) -> List[Dict]:
        """
        Get all conditions associated with a specific gene.
        
        Args:
            gene: Gene symbol to search for
            
        Returns:
            List[Dict]: List of conditions associated with the gene
        """
        return [
            condition for condition in self.conditions_data.get("conditions", [])
            if condition.get("genetic_data", {}).get("gene") == gene
        ]
    
    def get_condition_by_variant_id(self, variant_id: str) -> List[Dict]:
        """
        Get all conditions associated with a specific variant ID (rsID).
        
        Args:
            variant_id: Variant ID (rsID) to search for
            
        Returns:
            List[Dict]: List of conditions associated with the variant ID
        """
        return [
            condition for condition in self.conditions_data.get("conditions", [])
            if condition.get("genetic_data", {}).get("variant_id") == variant_id
        ]
    
    def get_conditions_by_risk(self, risk: str) -> List[Dict]:
        """
        Get all conditions with a specific risk level.
        
        Args:
            risk: Risk level to search for
            
        Returns:
            List[Dict]: List of conditions with the specified risk level
        """
        return [
            condition for condition in self.conditions_data.get("conditions", [])
            if condition.get("status", {}).get("risk") == risk
        ]
    
    def get_pathogenic_conditions(self) -> List[Dict]:
        """
        Get all pathogenic or likely pathogenic conditions.
        
        Returns:
            List[Dict]: List of pathogenic or likely pathogenic conditions
        """
        return [
            condition for condition in self.conditions_data.get("conditions", [])
            if condition.get("status", {}).get("risk") in ["Pathogenic", "Likely Pathogenic"]
        ]
    
    def annotate_variants(self, variants_df):
        """
        Annotate a variants DataFrame with clinical information.
        
        Args:
            variants_df: DataFrame of variants
            
        Returns:
            DataFrame: Annotated variants DataFrame
        """
        if variants_df is None or variants_df.empty:
            return variants_df
        
        # Create a copy to avoid modifying the original
        annotated_df = variants_df.copy()
        
        # Add clinical annotation columns
        annotated_df['clinical_significance'] = None
        annotated_df['condition_name'] = None
        annotated_df['condition_description'] = None
        annotated_df['rcv'] = None
        
        # Annotate each variant
        for idx, row in annotated_df.iterrows():
            gene = row.get('gene')
            rsid = row.get('rsid')
            
            # Look for matching conditions
            matching_conditions = []
            if gene and rsid:
                # Try exact match on variant_id first
                matching_conditions = self.get_condition_by_variant_id(rsid)
                
                # If no match, try gene match
                if not matching_conditions:
                    matching_conditions = self.get_condition_by_gene(gene)
            
            # Apply annotations if we found matches
            if matching_conditions:
                # Use the first match for now (could be enhanced to handle multiple matches)
                condition = matching_conditions[0]
                annotated_df.at[idx, 'clinical_significance'] = condition.get('status', {}).get('risk')
                annotated_df.at[idx, 'condition_name'] = condition.get('name')
                annotated_df.at[idx, 'condition_description'] = condition.get('description')
                annotated_df.at[idx, 'rcv'] = condition.get('genetic_data', {}).get('rcv')
        
        return annotated_df
    
    def generate_clinical_report(self, variants_df, output_dir: Union[str, Path]) -> str:
        """
        Generate a clinical report based on the variants and conditions data.
        
        Args:
            variants_df: DataFrame of variants
            output_dir: Directory to save the report
            
        Returns:
            str: Path to the generated report
        """
        if variants_df is None or variants_df.empty:
            log.warning("No variants provided for clinical report generation")
            return None
        
        # Annotate variants with clinical information
        annotated_df = self.annotate_variants(variants_df)
        
        # Prepare report content
        report_content = "# Clinical Variant Report\n\n"
        
        # Add summary section
        report_content += "## Summary\n\n"
        
        pathogenic_count = annotated_df[annotated_df['clinical_significance'].isin(
            ['Pathogenic', 'Likely Pathogenic'])].shape[0]
        
        benign_count = annotated_df[annotated_df['clinical_significance'].isin(
            ['Benign', 'Likely Benign'])].shape[0]
        
        vus_count = annotated_df[annotated_df['clinical_significance'] == 'Uncertain Significance'].shape[0]
        
        report_content += f"- **Pathogenic/Likely Pathogenic Variants**: {pathogenic_count}\n"
        report_content += f"- **Benign/Likely Benign Variants**: {benign_count}\n"
        report_content += f"- **Variants of Uncertain Significance**: {vus_count}\n\n"
        
        # Add detailed findings section
        report_content += "## Detailed Findings\n\n"
        
        # Group by clinical significance
        if pathogenic_count > 0:
            report_content += "### Pathogenic/Likely Pathogenic Variants\n\n"
            pathogenic_variants = annotated_df[annotated_df['clinical_significance'].isin(
                ['Pathogenic', 'Likely Pathogenic'])]
            
            for _, row in pathogenic_variants.iterrows():
                report_content += f"#### {row.get('gene')} - {row.get('rsid')}\n\n"
                report_content += f"- **Condition**: {row.get('condition_name')}\n"
                report_content += f"- **Description**: {row.get('condition_description')}\n"
                report_content += f"- **Genotype**: {row.get('genotype')}\n"
                report_content += f"- **Clinical Significance**: {row.get('clinical_significance')}\n"
                report_content += f"- **ClinVar ID**: {row.get('rcv')}\n\n"
        
        # Save the report
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        report_path = output_dir / "clinical_report.md"
        with open(report_path, 'w') as f:
            f.write(report_content)
        
        log.info(f"Generated clinical report: {report_path}")
        return str(report_path)


def create_example_conditions() -> Dict:
    """
    Create example conditions data for testing.
    
    Returns:
        Dict: Example conditions data
    """
    return {
        "conditions": [
            {
                "name": "CHRM2-Related Cognitive Function",
                "description": "Variants in CHRM2 associated with executive function, memory, and attention",
                "symptoms": [
                    "Reduced executive function",
                    "Memory impairment",
                    "Attention deficits"
                ],
                "status": {
                    "risk": "Risk Factor",
                    "confidence": "Medium",
                    "classification": "P"
                },
                "genetic_data": {
                    "gene": "CHRM2",
                    "variant_id": "rs2350780",
                    "rcv": "RCV000123456",
                    "genotype": "G/A"
                },
                "risk_assessment": {
                    "frequency": 0.15,
                    "version": "GRCh38"
                }
            },
            {
                "name": "DRD2-Related Dopamine Modulation",
                "description": "Variants in DRD2 associated with dopamine modulation and cognitive flexibility",
                "symptoms": [
                    "Reduced cognitive flexibility",
                    "Altered dopamine signaling"
                ],
                "status": {
                    "risk": "Risk Factor",
                    "confidence": "Medium",
                    "classification": "P"
                },
                "genetic_data": {
                    "gene": "DRD2",
                    "variant_id": "rs6277",
                    "rcv": "RCV000789012",
                    "genotype": "G/A"
                },
                "risk_assessment": {
                    "frequency": 0.25,
                    "version": "GRCh38"
                }
            },
            {
                "name": "Usher Syndrome Type II",
                "description": "A genetic disorder characterized by hearing loss and progressive vision loss",
                "symptoms": [
                    "Hearing loss",
                    "Progressive vision loss",
                    "Night blindness",
                    "Loss of peripheral vision"
                ],
                "status": {
                    "risk": "Pathogenic",
                    "confidence": "High",
                    "classification": "D"
                },
                "genetic_data": {
                    "gene": "ADGRV1",
                    "variant_id": "rs575602255",
                    "rcv": "RCV000345678",
                    "genotype": "A/G"
                },
                "risk_assessment": {
                    "frequency": 0.001,
                    "version": "GRCh38"
                }
            },
            {
                "name": "NBIA-4 Neurodegeneration",
                "description": "Neurodegeneration with Brain Iron Accumulation type 4, characterized by motor and cognitive decline",
                "symptoms": [
                    "Progressive motor dysfunction",
                    "Cognitive decline",
                    "Spasticity",
                    "Dystonia"
                ],
                "status": {
                    "risk": "Pathogenic",
                    "confidence": "High",
                    "classification": "D"
                },
                "genetic_data": {
                    "gene": "C19orf12",
                    "variant_id": "rs146170087",
                    "rcv": "RCV000567890",
                    "genotype": "T/C"
                },
                "risk_assessment": {
                    "frequency": 0.0005,
                    "version": "GRCh38"
                }
            }
        ]
    }


if __name__ == "__main__":
    # Example usage
    logging.basicConfig(level=logging.INFO)
    
    # Create a clinical variant manager
    manager = ClinicalVariantManager()
    
    # Create example data
    example_data = create_example_conditions()
    
    # Save example data to file
    output_dir = Path("examples")
    output_dir.mkdir(exist_ok=True)
    
    with open(output_dir / "example_conditions.json", 'w') as f:
        json.dump(example_data, f, indent=2)
    
    print(f"Example conditions saved to {output_dir / 'example_conditions.json'}")
