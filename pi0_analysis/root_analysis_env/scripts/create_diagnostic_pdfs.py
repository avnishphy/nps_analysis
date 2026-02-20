#!/usr/bin/env python3
"""
Create a single PDF file from ROOT diagnostic files.
Extracts all plots from all diagnostics_run*.root files and arranges them 3x3 per page.
"""

import os
import sys
import glob
import subprocess
from pathlib import Path
from typing import List, Tuple
import tempfile
import shutil

try:
    import ROOT
    from ROOT import TFile, TCanvas
except ImportError:
    print("ERROR: ROOT Python bindings not available. Please set up ROOT environment.")
    sys.exit(1)

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    print("ERROR: PIL/Pillow not available. Install with: pip install Pillow")
    sys.exit(1)

try:
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfgen import canvas as pdf_canvas
    from reportlab.lib.units import inch
except ImportError:
    print("ERROR: reportlab not available. Install with: pip install reportlab")
    sys.exit(1)


class DiagnosticPDFGenerator:
    """Generate single PDF from ROOT diagnostic files."""
    
    def __init__(self, output_dir: str = "output/plots/x60_4b"):
        self.output_dir = Path(output_dir).resolve()  # Resolve to absolute path
        self.temp_dir = tempfile.mkdtemp(prefix="root_diagnostics_")
        self.page_size = letter  # 8.5 x 11 inches
        self.margin = 0.25 * inch
        self.plots_per_page = 9  # 3x3 grid
        self.grid_cols = 3
        self.grid_rows = 3
        
        print(f"Output directory: {self.output_dir}")
        print(f"Temp directory: {self.temp_dir}")
    
    def find_diagnostic_files(self) -> List[Tuple[int, Path]]:
        """Find all diagnostics_run*.root files and return sorted list."""
        pattern = str(self.output_dir / "diagnostics_run*.root")
        files = glob.glob(pattern)
        
        # Extract run numbers and sort
        run_files = []
        for f in files:
            try:
                run_num = int(Path(f).stem.split('_run')[1])
                run_files.append((run_num, Path(f)))
            except (IndexError, ValueError):
                print(f"Warning: Could not parse run number from {f}")
        
        run_files.sort(key=lambda x: x[0])
        return run_files
    
    def extract_canvases_from_root(self, root_file: Path, run_num: int) -> List[Tuple[str, Path, int]]:
        """Extract all canvases from ROOT file and save as PNG images.
        Returns list of (canvas_name, png_path, run_num) tuples."""
        canvas_images = []
        
        try:
            # Disable ROOT GUI
            ROOT.gROOT.SetBatch(True)
            
            # Open ROOT file
            tfile = TFile(str(root_file), "READ")
            if not tfile or tfile.IsZombie():
                print(f"ERROR: Cannot open {root_file}")
                return canvas_images
            
            # Get list of all objects in file
            key_list = tfile.GetListOfKeys()
            
            for i in range(key_list.GetSize()):
                key = key_list.At(i)
                obj_name = key.GetName()
                obj_class = key.GetClassName()
                
                # Look for TCanvas objects
                if obj_class == "TCanvas":
                    obj = tfile.Get(obj_name)
                    if obj:
                        # Save canvas as PNG
                        png_path = Path(self.temp_dir) / f"run{run_num}_{obj_name}.png"
                        obj.SaveAs(str(png_path))
                        
                        if png_path.exists():
                            canvas_images.append((obj_name, png_path, run_num))
                            print(f"  Extracted: {obj_name}")
            
            tfile.Close()
        
        except Exception as e:
            print(f"ERROR processing {root_file}: {e}")
        
        return canvas_images
    
    def create_grid_page(self, images: List[Tuple[str, Path, int]], page_num: int) -> Path:
        """Create a PDF page with 3x3 grid of images with run and canvas info."""
        try:
            # Page dimensions in pixels (at 150 DPI for good quality)
            dpi = 150
            page_width_px = int(self.page_size[0] / inch * dpi)
            page_height_px = int(self.page_size[1] / inch * dpi)
            
            # Calculate plot area (accounting for margins)
            margin_px = int(self.margin / inch * dpi)
            plot_width = (page_width_px - 2 * margin_px) // self.grid_cols
            plot_height = (page_height_px - 4 * margin_px - 60) // self.grid_rows  # 60px for title
            
            # Create blank page
            page = Image.new('RGB', (page_width_px, page_height_px), 'white')
            draw = ImageDraw.Draw(page)
            
            # Try to load a font (fallback to default if not available)
            try:
                title_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 20)
                label_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 10)
            except:
                title_font = ImageFont.load_default()
                label_font = ImageFont.load_default()
            
            # Draw title
            title = f"Diagnostic Plots - Page {page_num}"
            title_bbox = draw.textbbox((0, 0), title, font=title_font)
            title_width = title_bbox[2] - title_bbox[0]
            title_x = (page_width_px - title_width) // 2
            draw.text((title_x, margin_px // 2), title, fill='black', font=title_font)
            
            # Place images in grid
            for idx, (canvas_name, img_path, run_num) in enumerate(images):
                row = idx // self.grid_cols
                col = idx % self.grid_cols
                
                # Calculate position
                x = margin_px + col * plot_width
                y = margin_px * 2 + 60 + row * plot_height
                
                # Load and resize image
                try:
                    img = Image.open(img_path)
                    img.thumbnail((plot_width - 10, plot_height - 35), Image.Resampling.LANCZOS)
                    
                    # Center image within grid cell
                    img_x = x + (plot_width - img.width) // 2
                    img_y = y + (plot_height - img.height) // 2 - 10
                    
                    # Paste image
                    if img.mode == 'RGBA':
                        page.paste(img, (int(img_x), int(img_y)), img)
                    else:
                        page.paste(img, (int(img_x), int(img_y)))
                    
                    # Draw border around cell
                    draw.rectangle(
                        [(x, y), (x + plot_width, y + plot_height)],
                        outline='gray',
                        width=1
                    )
                    
                    # Draw label with run number and canvas name
                    label = f"Run {run_num}: {canvas_name}"
                    draw.text((x + 5, y + plot_height - 25), label, fill='darkblue', font=label_font)
                
                except Exception as e:
                    print(f"WARNING: Could not process image {img_path}: {e}")
            
            # Save page as temporary image
            temp_img_path = Path(self.temp_dir) / f"page_{page_num}.png"
            page.save(str(temp_img_path), 'PNG')
            
            return temp_img_path
        
        except Exception as e:
            print(f"ERROR creating grid page: {e}")
            return None
    
    def convert_images_to_pdf(self, images: List[Path], output_pdf: Path) -> bool:
        """Convert image pages to PDF using ImageMagick."""
        try:
            if not images:
                print(f"No images to convert")
                return False
            
            # Use ImageMagick to convert images to PDF
            image_files = [str(img) for img in images]
            cmd = ['convert'] + image_files + ['-quality', '85', str(output_pdf)]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"WARNING: ImageMagick conversion failed: {result.stderr}")
                # Try alternative method with reportlab
                return self._convert_with_reportlab(images, output_pdf)
            
            return True
        
        except Exception as e:
            print(f"ERROR converting images to PDF: {e}")
            return self._convert_with_reportlab(images, output_pdf)
    
    def _convert_with_reportlab(self, images: List[Path], output_pdf: Path) -> bool:
        """Fallback PDF creation using reportlab."""
        try:
            c = pdf_canvas.Canvas(str(output_pdf), pagesize=self.page_size)
            
            for img_path in images:
                c.drawImage(str(img_path), 0, 0, 
                           width=self.page_size[0], 
                           height=self.page_size[1])
                c.showPage()
            
            c.save()
            return True
        
        except Exception as e:
            print(f"ERROR with reportlab conversion: {e}")
            return False
    
    def run(self) -> int:
        """Main execution."""
        print("\n" + "="*70)
        print("ROOT Diagnostic PDF Generator (Single Combined PDF)")
        print("="*70)
        print(f"Looking for ROOT files in: {self.output_dir}")
        
        # Find all diagnostic files
        run_files = self.find_diagnostic_files()
        
        if not run_files:
            print(f"ERROR: No diagnostic files found in {self.output_dir}")
            print(f"Try running with absolute path or from the workspace root directory")
            self.cleanup()
            return 1
        
        print(f"\nFound {len(run_files)} diagnostic files:")
        for run_num, fpath in run_files:
            print(f"  Run {run_num}: {fpath.name}")
        
        # Extract all canvases from all runs
        print("\nExtracting canvases from all runs...")
        all_canvases = []
        for run_num, root_file in run_files:
            print(f"\nProcessing Run {run_num}...")
            canvases = self.extract_canvases_from_root(root_file, run_num)
            all_canvases.extend(canvases)
            print(f"  Found {len(canvases)} canvases")
        
        if not all_canvases:
            print(f"ERROR: No canvases found in any diagnostic files")
            self.cleanup()
            return 1
        
        print(f"\nTotal canvases extracted: {len(all_canvases)}")
        
        # Sort by run number, then by canvas name
        all_canvases.sort(key=lambda x: (x[2], x[0]))
        
        # Create grid pages
        print("\nCreating PDF pages...")
        page_images = []
        page_num = 1
        
        for i, (canvas_name, img_path, run_num) in enumerate(all_canvases):
            page_images.append((canvas_name, img_path, run_num))
            
            # Check if we've filled a page or reached the end
            if len(page_images) == self.plots_per_page or i == len(all_canvases) - 1:
                print(f"Creating page {page_num} with {len(page_images)} images...")
                page_path = self.create_grid_page(page_images, page_num)
                if page_path:
                    page_num += 1
                
                page_images = []
        
        # Create PDF from pages
        page_images = sorted(
            Path(self.temp_dir).glob("page_*.png"),
            key=lambda x: int(x.stem.split('_')[1])
        )
        
        output_pdf = self.output_dir / "all_diagnostics.pdf"
        
        if self.convert_images_to_pdf(page_images, output_pdf):
            print(f"\n✓ Successfully created {output_pdf.name}")
            print(f"  Location: {output_pdf}")
            print(f"  Total pages: {len(page_images)}")
            return 0
        else:
            print(f"\n✗ Failed to create PDF")
            return 1
    
    def cleanup(self):
        """Remove temporary directory."""
        try:
            shutil.rmtree(self.temp_dir)
            print(f"\nCleaned up temporary directory")
        except Exception as e:
            print(f"Warning: Could not remove temp directory {self.temp_dir}: {e}")


def main():
    """Entry point."""
    # Get output directory from command line or use default
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "output/plots/x60_4b"
    
    generator = DiagnosticPDFGenerator(output_dir)
    result = generator.run()
    generator.cleanup()
    return result


if __name__ == "__main__":
    sys.exit(main())

    
    def find_diagnostic_files(self) -> List[Tuple[int, Path]]:
        """Find all diagnostics_run*.root files and return sorted list."""
        pattern = str(self.output_dir / "diagnostics_run*.root")
        files = glob.glob(pattern)
        
        # Extract run numbers and sort
        run_files = []
        for f in files:
            try:
                run_num = int(Path(f).stem.split('_run')[1])
                run_files.append((run_num, Path(f)))
            except (IndexError, ValueError):
                print(f"Warning: Could not parse run number from {f}")
        
        run_files.sort(key=lambda x: x[0])
        return run_files
    
    def extract_canvases_from_root(self, root_file: Path) -> List[Tuple[str, Path]]:
        """Extract all canvases from ROOT file and save as PNG images."""
        canvas_images = []
        
        try:
            # Disable ROOT GUI
            ROOT.gROOT.SetBatch(True)
            
            # Open ROOT file
            tfile = TFile(str(root_file), "READ")
            if not tfile or tfile.IsZombie():
                print(f"ERROR: Cannot open {root_file}")
                return canvas_images
            
            # Get list of all objects in file
            key_list = tfile.GetListOfKeys()
            
            for i in range(key_list.GetSize()):
                key = key_list.At(i)
                obj_name = key.GetName()
                obj_class = key.GetClassName()
                
                # Look for TCanvas objects
                if obj_class == "TCanvas":
                    obj = tfile.Get(obj_name)
                    if obj:
                        # Save canvas as PNG
                        png_path = Path(self.temp_dir) / f"{obj_name}.png"
                        obj.SaveAs(str(png_path))
                        
                        if png_path.exists():
                            canvas_images.append((obj_name, png_path))
                            print(f"  Extracted: {obj_name}")
            
            tfile.Close()
        
        except Exception as e:
            print(f"ERROR processing {root_file}: {e}")
        
        return canvas_images
    
    def create_grid_page(self, images: List[Path], run_num: int, 
                        page_num: int, output_pdf: Path) -> bool:
        """Create a PDF page with 3x3 grid of images."""
        try:
            # Determine actual grid size (might be less than 3x3 on last page)
            num_images = len(images)
            actual_rows = (num_images + self.grid_cols - 1) // self.grid_cols
            
            # Page dimensions in pixels (at 150 DPI for good quality)
            dpi = 150
            page_width_px = int(self.page_size[0] / inch * dpi)
            page_height_px = int(self.page_size[1] / inch * dpi)
            
            # Calculate plot area (accounting for margins)
            margin_px = int(self.margin / inch * dpi)
            plot_width = (page_width_px - 2 * margin_px) // self.grid_cols
            plot_height = (page_height_px - 4 * margin_px - 60) // self.grid_rows  # 60px for title
            
            # Create blank page
            page = Image.new('RGB', (page_width_px, page_height_px), 'white')
            draw = ImageDraw.Draw(page)
            
            # Try to load a font (fallback to default if not available)
            try:
                title_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 20)
                label_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12)
            except:
                title_font = ImageFont.load_default()
                label_font = ImageFont.load_default()
            
            # Draw title
            title = f"Run {run_num} - Diagnostics (Page {page_num})"
            title_bbox = draw.textbbox((0, 0), title, font=title_font)
            title_width = title_bbox[2] - title_bbox[0]
            title_x = (page_width_px - title_width) // 2
            draw.text((title_x, margin_px // 2), title, fill='black', font=title_font)
            
            # Place images in grid
            for idx, img_path in enumerate(images):
                row = idx // self.grid_cols
                col = idx % self.grid_cols
                
                # Calculate position
                x = margin_px + col * plot_width
                y = margin_px * 2 + 60 + row * plot_height
                
                # Load and resize image
                try:
                    img = Image.open(img_path)
                    img.thumbnail((plot_width - 10, plot_height - 30), Image.Resampling.LANCZOS)
                    
                    # Center image within grid cell
                    img_x = x + (plot_width - img.width) // 2
                    img_y = y + (plot_height - img.height) // 2
                    
                    # Paste image
                    if img.mode == 'RGBA':
                        page.paste(img, (int(img_x), int(img_y)), img)
                    else:
                        page.paste(img, (int(img_x), int(img_y)))
                    
                    # Draw border around cell
                    draw.rectangle(
                        [(x, y), (x + plot_width, y + plot_height)],
                        outline='gray',
                        width=1
                    )
                    
                    # Draw filename label
                    label = Path(img_path).stem
                    draw.text((x + 5, y + plot_height - 25), label, fill='gray', font=label_font)
                
                except Exception as e:
                    print(f"WARNING: Could not process image {img_path}: {e}")
            
            # Save page as temporary image
            temp_img_path = Path(self.temp_dir) / f"page_run{run_num}_{page_num}.png"
            page.save(str(temp_img_path), 'PNG')
            
            return True
        
        except Exception as e:
            print(f"ERROR creating grid page: {e}")
            return False
    
    def convert_images_to_pdf(self, images: List[Path], run_num: int, 
                             output_pdf: Path) -> bool:
        """Convert image pages to PDF using ImageMagick."""
        try:
            if not images:
                print(f"No images to convert for run {run_num}")
                return False
            
            # Use ImageMagick to convert images to PDF
            image_files = [str(img) for img in images]
            cmd = ['convert'] + image_files + ['-quality', '85', str(output_pdf)]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"WARNING: ImageMagick conversion failed: {result.stderr}")
                # Try alternative method with reportlab
                return self._convert_with_reportlab(images, run_num, output_pdf)
            
            return True
        
        except Exception as e:
            print(f"ERROR converting images to PDF: {e}")
            return self._convert_with_reportlab(images, run_num, output_pdf)
    
    def _convert_with_reportlab(self, images: List[Path], run_num: int,
                               output_pdf: Path) -> bool:
        """Fallback PDF creation using reportlab."""
        try:
            c = pdf_canvas.Canvas(str(output_pdf), pagesize=self.page_size)
            
            for img_path in images:
                c.drawImage(str(img_path), 0, 0, 
                           width=self.page_size[0], 
                           height=self.page_size[1])
                c.showPage()
            
            c.save()
            return True
        
        except Exception as e:
            print(f"ERROR with reportlab conversion: {e}")
            return False
    
    def process_run(self, run_num: int, root_file: Path) -> bool:
        """Process a single run and create PDF."""
        print(f"\nProcessing Run {run_num}...")
        
        # Extract canvases
        canvases = self.extract_canvases_from_root(root_file)
        
        if not canvases:
            print(f"No canvases found in {root_file}")
            return False
        
        print(f"Found {len(canvases)} canvases")
        
        # Sort canvases by name for consistent ordering
        canvases.sort(key=lambda x: x[0])
        
        # Create grid pages
        page_images = []
        page_num = 1
        
        for i, (name, img_path) in enumerate(canvases):
            page_images.append(img_path)
            
            # Check if we've filled a page or reached the end
            if len(page_images) == self.plots_per_page or i == len(canvases) - 1:
                print(f"Creating page {page_num} with {len(page_images)} images...")
                page_path = Path(self.temp_dir) / f"page_run{run_num}_{page_num}.png"
                
                if self.create_grid_page(page_images, run_num, page_num, page_path):
                    page_num += 1
                
                page_images = []
        
        # Create PDF from pages
        page_images = sorted(
            Path(self.temp_dir).glob(f"page_run{run_num}_*.png"),
            key=lambda x: int(x.stem.split('_')[-1])
        )
        
        output_pdf = self.output_dir / f"diagnostics_run{run_num}.pdf"
        
        if self.convert_images_to_pdf(page_images, run_num, output_pdf):
            print(f"✓ Successfully created {output_pdf.name}")
            return True
        else:
            print(f"✗ Failed to create PDF for run {run_num}")
            return False
    
    def cleanup(self):
        """Remove temporary directory."""
        try:
            shutil.rmtree(self.temp_dir)
            print(f"\nCleaned up temporary directory")
        except Exception as e:
            print(f"Warning: Could not remove temp directory {self.temp_dir}: {e}")
    
    def run(self) -> int:
        """Main execution."""
        print("\n" + "="*70)
        print("ROOT Diagnostic PDF Generator")
        print("="*70)
        print(f"Looking for ROOT files in: {self.output_dir}")
        
        # Find all diagnostic files
        run_files = self.find_diagnostic_files()
        
        if not run_files:
            print(f"ERROR: No diagnostic files found in {self.output_dir}")
            print(f"Try running with absolute path or from the workspace root directory")
            self.cleanup()
            return 1
        
        print(f"\nFound {len(run_files)} diagnostic files:")
        for run_num, fpath in run_files:
            print(f"  Run {run_num}: {fpath.name}")
        
        # Process each run
        successful = 0
        failed = 0
        
        for run_num, root_file in run_files:
            if self.process_run(run_num, root_file):
                successful += 1
            else:
                failed += 1
        
        # Summary
        print("\n" + "="*70)
        print(f"Summary: {successful} successful, {failed} failed")
        print("="*70 + "\n")
        
        # Cleanup
        self.cleanup()
        
        return 0 if failed == 0 else 1


def main():
    """Entry point."""
    # Get output directory from command line or use default
    output_dir = sys.argv[1] if len(sys.argv) > 1 else "output/plots/x60_4b"
    
    generator = DiagnosticPDFGenerator(output_dir)
    return generator.run()


if __name__ == "__main__":
    sys.exit(main())
