"""
Module for analyzing Structural Variants (SVs), specifically Copy Number Variations (CNVs),
within the RCCX locus (chr6:p21.3).

Focuses on identifying potential CNVs reported within VCF files.
"""

import logging
from typing import Tuple, Dict, Any, List, Optional

# Using GRCh38/hg38 coordinates for the broader RCCX region
# Spanning approximately from STK19/RP1 to ZKSCAN16/DOM3Z including C4A/B, CYP21A2, TNXB
# chr6:31,950,000-32,150,000 gives some buffer around the core ~190kb module.
RCCX_REGION_GRCH38: Dict[str, Any] = {
    "chrom": "chr6",
    "start": 31_950_000,
    "end": 32_150_000,
    "name": "RCCX Locus (GRCh38)"
}

# TODO: Add support for GRCh37 coordinates if necessary

logger = logging.getLogger(__name__)

def get_rccx_region(reference: str = "GRCh38") -> Optional[Dict[str, Any]]:
    """Returns the defined genomic coordinates for the RCCX region."""
    if reference.upper() in ["GRCH38", "HG38"]:
        # Ensure chromosome is in the expected format (e.g., '6' vs 'chr6')
        # This might need adjustment based on input VCF contig format
        region = RCCX_REGION_GRCH38.copy()
        # Example adjustment (needs confirmation based on VCF headers):
        # if vcf_uses_numeric_chrom:
        #    region["chrom"] = "6"
        return region
    # elif reference.upper() in ["GRCH37", "HG19"]:
        # return RCCX_REGION_GRCH37.copy() # Add hg19 coords if needed
    else:
        logger.error(f"Unsupported reference genome: {reference}")
        return None

# Placeholder for the main analysis function
def analyze_rccx_cnv_from_vcf(
    vcf_path: str,
    reference: str = "GRCh38"
) -> List[Dict[str, Any]]:
    """
    Analyzes a VCF file for potential CNVs within the RCCX region.

    Args:
        vcf_path: Path to the input VCF file (must be indexed).
        reference: Reference genome build (e.g., "GRCh38").

    Returns:
        A list of dictionaries, each representing a potential RCCX CNV found.
        Returns an empty list if no relevant CNVs are found or errors occur.
    """
    results: List[Dict[str, Any]] = []
    rccx_region = get_rccx_region(reference)
    if not rccx_region:
        return results # Error logged in get_rccx_region

    try:
        import pysam # Optional dependency, ensure it's managed by Poetry if used heavily
    except ImportError:
        logger.error("pysam library is required for RCCX CNV analysis from VCF. pip install pysam")
        # Or add to pyproject.toml: poetry add pysam
        return results

    try:
        vcf_file = pysam.VariantFile(vcf_path)

        # Check VCF contig format matches expected region format
        vcf_contigs = list(vcf_file.header.contigs)
        target_chrom = rccx_region["chrom"]
        if target_chrom not in vcf_contigs:
            # Attempt conversion (e.g., chr6 -> 6 or 6 -> chr6)
            if target_chrom.startswith("chr") and target_chrom[3:] in vcf_contigs:
                target_chrom = target_chrom[3:]
                logger.warning(f"Adjusting target chromosome to '{target_chrom}' based on VCF header.")
            elif not target_chrom.startswith("chr") and f"chr{target_chrom}" in vcf_contigs:
                target_chrom = f"chr{target_chrom}"
                logger.warning(f"Adjusting target chromosome to '{target_chrom}' based on VCF header.")
            else:
                 logger.error(f"Target chromosome '{rccx_region['chrom']}' not found in VCF contigs: {vcf_contigs}. Cannot perform RCCX analysis.")
                 vcf_file.close()
                 return results
        
        # Fetch variants within the RCCX region
        fetched_variants = vcf_file.fetch(
            target_chrom,
            rccx_region["start"],
            rccx_region["end"]
        )

        for record in fetched_variants:
            sv_info = parse_sv_from_vcf_record(record)
            if sv_info:
                # Basic overlap check (can be refined)
                # Ensure the reported SV overlaps the core RCCX region
                if sv_info['chrom'] == target_chrom and \
                   max(sv_info['start'], rccx_region['start']) < min(sv_info['end'], rccx_region['end']):
                    results.append(sv_info)

        vcf_file.close()

    except FileNotFoundError:
        logger.error(f"VCF file not found: {vcf_path}")
    except ValueError as e:
        logger.error(f"Error processing VCF {vcf_path}. Is it indexed (e.g., with tabix)? Error: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred during RCCX CNV analysis: {e}", exc_info=True)

    if results:
         logger.info(f"Found {len(results)} potential SV/CNV records within the RCCX region in {vcf_path}")
    else:
         logger.info(f"No potential SV/CNV records identified in the RCCX region within {vcf_path}")

    return results


def parse_sv_from_vcf_record(record: 'pysam.VariantRecord') -> Optional[Dict[str, Any]]:
    """
    Attempts to parse Structural Variant information from a VCF record.

    Checks for standard INFO fields (SVTYPE, END, SVLEN) and ALT allele formats.

    Args:
        record: A pysam VariantRecord object.

    Returns:
        A dictionary containing SV information (type, chrom, start, end, length, raw_record)
        if found, otherwise None.
    """
    info = record.info
    alt_alleles = record.alts or []

    sv_type = info.get('SVTYPE')
    sv_end = info.get('END')
    sv_len = info.get('SVLEN') # Often negative for deletions

    # Try inferring from ALT alleles if INFO fields are missing
    if not sv_type:
        for alt in alt_alleles:
            if alt in ["<DEL>", "<CN0>"]:
                sv_type = "DEL"
                break
            elif alt in ["<DUP>", "<DUP:TANDEM>", "<CN2+>"]: # Simplification
                sv_type = "DUP"
                break
            elif alt == "<CNV>":
                sv_type = "CNV" # Generic CNV
                break
            # Could add more complex ALT parsing if needed

    # If we identified an SV type, try to get coordinates
    if sv_type:
        start = record.pos # VCF POS is 1-based start
        chrom = record.chrom

        # VCF spec: END is required for SVs unless ALT defines length implicitly (e.g. symbolic)
        # SVLEN is also common. Prioritize END if available.
        end = sv_end
        if end is None:
            # Try to calculate end from SVLEN if present
            if sv_len is not None:
                # SVLEN is difference: length = ALT - REF. Negative for DEL.
                # For simple DEL/DUP, end = start + abs(sv_len) ? No, this isn't quite right.
                # VCF POS is the base *before* the event for DEL/INS.
                # For simple DEL: end = start + abs(svlen) -1 ? No...
                # For simple DUP: end = start + svlen -1 ? No...
                # Let's require END for now for reliable coordinates unless it's simple <DEL>/<DUP>
                # or point to the record itself if END missing.
                # A simple DEL's length is len(REF) - len(ALT) = len(REF) - 1
                # A simple DUP (INS) length is len(ALT) - len(REF) = len(ALT) - 1
                 if sv_type == "DEL" and sv_len is not None:
                     # If SVLEN is negative, its absolute value is the length of the deletion
                     end = start + abs(sv_len) # Approx end if SVLEN is length of deleted seq
                 # DUPs are more complex to get precise end without END field.

        # Fallback if END still not determined
        if end is None:
             # For symbolic alleles like <DEL>, often the END tag IS the end coordinate.
             # If END is missing, the precise boundaries are uncertain from this record alone.
             logger.debug(f"SV record at {chrom}:{start} has type {sv_type} but missing END INFO field. Boundaries uncertain.")
             # We could potentially make assumptions based on REF/ALT lengths for non-symbolic SVs,
             # but symbolic ones like <DEL>/<DUP> absolutely need END.
             # Let's skip records with uncertain boundaries for now.
             return None


        # Ensure end is valid
        if not isinstance(end, int) or end < start:
             logger.warning(f"Invalid or missing END coordinate ({end}) for SV at {chrom}:{start}. Skipping.")
             return None

        return {
            "type": sv_type,
            "chrom": chrom,
            "start": start, # VCF 1-based
            "end": end,     # VCF 1-based
            "length": abs(sv_len) if sv_len is not None else (end - start + 1), # Approx length
            "record_details": str(record).strip() # Provide raw record context
        }

    return None

# Example Usage (if run directly)
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    # Replace with a path to a test VCF file containing SVs in the RCCX region
    # Ensure the VCF is indexed (e.g., test.vcf.gz and test.vcf.gz.tbi)
    test_vcf = "path/to/your/test_sv.vcf.gz"
    
    logger.info(f"Analyzing VCF: {test_vcf}")
    found_cnvs = analyze_rccx_cnv_from_vcf(test_vcf)

    if found_cnvs:
        logger.info("\n--- Found RCCX CNVs/SVs ---")
        for cnv in found_cnvs:
            logger.info(f"  Type: {cnv['type']}, Region: {cnv['chrom']}:{cnv['start']}-{cnv['end']}, Approx Length: {cnv['length']}")
            # logger.info(f"    Record: {cnv['record_details']}") # Uncomment for full record details
    else:
        logger.info("No relevant CNVs/SVs found in the specified VCF within the RCCX region.")

    # Example of getting region
    region = get_rccx_region()
    if region:
        logger.info(f"\nUsing RCCX region: {region['name']} {region['chrom']}:{region['start']}-{region['end']}")
