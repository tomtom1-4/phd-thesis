import sys
import fitz  # PyMuPDF
from PIL import Image
import io

def crop_pdf(input_pdf_path, output_pdf_path):
    """
    Crops a PDF file to eliminate white spaces around the visible content.

    :param input_pdf_path: Path to the input PDF file.
    :param output_pdf_path: Path where the cropped PDF will be saved.
    """
    with fitz.open(input_pdf_path) as doc:
        for page_number, page in enumerate(doc, start=1):
            # Render the page to an image
            pix = page.get_pixmap(alpha=False)  # Render the page to an image
            img_data = pix.tobytes("ppm")  # Get the image data in PPM format

            # Open the image with Pillow
            img = Image.open(io.BytesIO(img_data))

            # Convert image to grayscale
            gray = img.convert("L")

            # Create a binary image (black and white)
            # Set all pixels equal to 255 (white) to 0, and all other pixels to 1
            bw = gray.point(lambda x: 0 if x == 255 else 1, '1')

            # Get the bounding box of non-zero regions in the binary image
            bbox = bw.getbbox()

            if not bbox:
                print(f"Page {page_number} has no visible content, skipping cropping.")
                continue

            # Calculate PDF coordinates from pixel coordinates
            # Image coordinates origin is at the top-left corner
            # PDF page coordinates origin is at the bottom-left corner
            img_width, img_height = bw.size
            x0_img, y0_img, x1_img, y1_img = bbox

            # Calculate the crop box in PDF coordinates
            page_width, page_height = page.rect.width, page.rect.height

            # Scale factors
            x_scale = page_width / img_width
            y_scale = page_height / img_height

            # Map image coordinates to PDF coordinates
            x0_pdf = x0_img * x_scale
            x1_pdf = x1_img * x_scale
            y0_pdf = y0_img * y_scale
            y1_pdf = y1_img * y_scale
            #y0_pdf = (img_height - y1_img) * y_scale
            #y1_pdf = (img_height - y0_img) * y_scale

            # Create the new crop box rectangle
            crop_rect = fitz.Rect(x0_pdf, y0_pdf, x1_pdf, y1_pdf)

            # Optionally, expand the crop rectangle by a small margin
            margin = 5  # Adjust as needed
            crop_rect.x0 = max(crop_rect.x0 - margin, 0)
            crop_rect.y0 = max(crop_rect.y0 - margin, 0)
            crop_rect.x1 = min(crop_rect.x1 + margin, page.rect.width)
            crop_rect.y1 = min(crop_rect.y1 + margin, page.rect.height)

            # Set the page's crop box
            page.set_cropbox(crop_rect)
            print(f"Cropped page {page_number} to {crop_rect}.")

        # Save the cropped PDF to a new file
        doc.save(output_pdf_path)
        print(f"Cropped PDF saved as '{output_pdf_path}'.")

if __name__ == "__main__":
    if (len(sys.argv) != 2 and len(sys.argv) != 3):
        print("Usage: python crop.py input.pdf")
        sys.exit(1)

    input_pdf = sys.argv[1]
    if (len(sys.argv)==2):
        output_pdf = f"cropped_{input_pdf}"
    elif (len(sys.argv)==3):
        output_pdf = sys.argv[2]

    crop_pdf(input_pdf, output_pdf)