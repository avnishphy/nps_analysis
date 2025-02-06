import ROOT
import os
import re
import glob
import gc
import psutil

def plot_all_runs_to_single_pdf(files, output_pdf):
    """
    Plot histograms from multiple ROOT files and save to a single PDF.
    
    Args:
        files (list): List of tuples containing (file_path, run_number).
        output_pdf (str): Path to the output PDF.
    """
    first_page = True  # Flag to handle the first page in the PDF

    for file_path, run_number in files:
        file1 = ROOT.TFile.Open(file_path)
        keys_file1 = file1.GetListOfKeys()

        for key1 in keys_file1:
            name1 = key1.GetName()

            # Check if key starts with "hdc" and is TH1F
            if name1.startswith("hdc") and name1.endswith("ddist"):
                obj1 = key1.ReadObj()

                if isinstance(obj1, ROOT.TH1F):
                    try:
                        # Create a canvas
                        canvas = ROOT.TCanvas("canvas", f"Run {run_number} - {name1}", 800, 600)
                        canvas.SetGrid()

                        # Set line colors and styles
                        obj1.SetLineColor(ROOT.kRed)

                        # Draw histogram
                        obj1.Draw()

                        # Add run number to the canvas
                        run_label = ROOT.TLatex()
                        run_label.SetNDC()
                        run_label.SetTextSize(0.03)
                        run_label.SetTextAlign(11)  # Align left
                        run_label.DrawLatex(0.1, 0.93, f"Run Number: {run_number}")

                        # Create a legend
                        legend = ROOT.TLegend(0.7, 0.2, 0.9, 0.1)
                        legend.SetBorderSize(0)
                        legend.AddEntry(obj1, "calibrated", "l")
                        legend.Draw()

                        # Save to PDF
                        if first_page:
                            canvas.Print(output_pdf + "(", "pdf")  # Start the PDF
                            first_page = False
                        else:
                            canvas.Print(output_pdf, "pdf")

                        print(f"Added {name1} from Run {run_number} to PDF")
                    except Exception as e:
                        print(f"Error processing histogram {name1} from Run {run_number}: {e}")
                    finally:
                        # Cleanup
                        canvas.Close()
                        del canvas
                        del legend
                        del obj1

        # Close the ROOT file
        file1.Close()

    # Close the PDF
    if not first_page:
        ROOT.TCanvas().Print(output_pdf + ")", "pdf")  # End the PDF
    print(f"All histograms from all runs saved to {output_pdf}")

def log_memory_usage():
    """Logs the current memory usage of the Python process."""
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    print(f"Memory usage: {memory_info.rss / (1024**2):.2f} MB")

def main():
    directory = "/cache/hallc/c-nps/analysis/pass1/replays/skim"
    file_pattern = os.path.join(directory, "nps_hms_skim_*")

    # Get all file paths matching the pattern
    filepaths = glob.glob(file_pattern)
    files_with_runs = []  # List to store (file_path, run_number)

    for filepath in filepaths:
        try:
            # Extract the run number from the filename
            runnumber = int(filepath.split("_")[3])  # Get the run number part (e.g., '6810')

            # Append the file path and run number to the list
            files_with_runs.append((filepath, runnumber))
        except (IndexError, ValueError):
            print(f"Skipping file with unexpected name format: {filepath}")
            continue

    # Sort the list by run number (ascending)
    files_with_runs.sort(key=lambda x: x[1])
    length = len(files_with_runs)
    print("No. of runs: ", length)

    # Output PDF file path
    output_pdf = "/u/group/nps/singhav/nps_analysis/det_calibrations/drift_chamber/dc_comparison_pass1_replay/dc_comparison_all_runs_1.pdf"

    # Generate the single PDF with all histograms
    try:
        plot_all_runs_to_single_pdf(files_with_runs, output_pdf)
    except Exception as e:
        print(f"Error generating PDF: {e}")

    # Log memory usage and force garbage collection
    log_memory_usage()
    gc.collect()

    print(f"All histograms saved to {output_pdf}")

if __name__ == "__main__":
    main()
