#!/bin/bash
# ============================================================================
# Create diagnostic PDFs from ROOT files
# Usage: ./scripts/create_diagnostic_pdfs.sh [output_dir]
# Example: ./scripts/create_diagnostic_pdfs.sh output/plots/x60_4b
# ============================================================================

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSPACE_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Use provided output_dir or default
OUTPUT_DIR="${1:-output/plots/x60_4b}"

# If relative path, make it relative to workspace root
if [[ ! "$OUTPUT_DIR" = /* ]]; then
    OUTPUT_DIR="$WORKSPACE_ROOT/$OUTPUT_DIR"
fi

echo "Workspace root: $WORKSPACE_ROOT"
echo "Creating diagnostic PDFs from ROOT files..."
echo "Output directory: $OUTPUT_DIR"
echo "Script: $SCRIPT_DIR"
echo ""

# Check if diagnostic files exist
if ! ls "$OUTPUT_DIR"/diagnostics_run*.root &>/dev/null; then
    echo "ERROR: No diagnostic files found in $OUTPUT_DIR"
    echo "Expected pattern: $OUTPUT_DIR/diagnostics_run*.root"
    exit 1
fi

echo "Found diagnostic files:"
ls -1 "$OUTPUT_DIR"/diagnostics_run*.root | head -5
echo "... ($(ls "$OUTPUT_DIR"/diagnostics_run*.root | wc -l) total)"
echo ""

# Check dependencies
python3 << 'PYEOF'
import sys
deps = ['PIL', 'reportlab']
missing = []

for dep in deps:
    try:
        if dep == 'PIL':
            from PIL import Image
        elif dep == 'reportlab':
            from reportlab.lib.pagesizes import letter
        print(f"✓ {dep} available")
    except ImportError:
        missing.append(dep)
        print(f"✗ {dep} NOT available")

if missing:
    print(f"\nWARNING: Missing optional dependencies: {', '.join(missing)}")
    print("Install with:")
    if 'PIL' in missing:
        print("  pip install Pillow")
    if 'reportlab' in missing:
        print("  pip install reportlab")
    print("(Script will still attempt to create PDFs using ImageMagick)")
else:
    print("\nAll optional dependencies satisfied!")
    sys.exit(0)
PYEOF

echo ""
echo "Starting PDF generation..."
echo "This may take a while..."
echo ""

# Run the Python script with absolute path
python3 "$SCRIPT_DIR/create_diagnostic_pdfs.py" "$OUTPUT_DIR"
exit $?
