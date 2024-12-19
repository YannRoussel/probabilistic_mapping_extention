import requests
import gzip
import shutil

def unzip_gz(input_file, output_file):
    """
    Unzip a .gz file.

    :param input_file: Path to the input .gz file.
    :param output_file: Path to save the decompressed output file.
    """
    with gzip.open(input_file, 'rb') as f_in:
        with open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

def download_file(url, filename):
    """
    Download a file from the given URL and save it with the specified filename.

    :param url: URL of the file to download.
    :param filename: Name of the file to save.
    :return: True if the download was successful, False otherwise.
    """
    try:
        # Send HTTP request to the URL
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors (e.g., 404, 500)

        # Open the file in binary write mode and write the contents
        with open(filename, 'wb') as file:
            file.write(response.content)

        return True  # Download successful
    except requests.RequestException as e:
        print(f"Error downloading file from {url}: {e}")
        return False  # Download failed

def main():
    """
    Main function to download processed single cell RNA counts.
    """
    url_list = [
        "https://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/20200513_Mouse_PatchSeq_Release_exon.v2.csv.tar",
        "https://data.nemoarchive.org/other/AIBS/AIBS_patchseq/transcriptome/scell/SMARTseq/processed/analysis/20200611/20200513_Mouse_PatchSeq_Release_intron.v2.csv.tar"
    ]
    filename_list = [
        "20200513_Mouse_PatchSeq_Release_exon.v2.csv.tar",
        "20200513_Mouse_PatchSeq_Release_intron.v2.csv.tar"
    ]

    for url, filename in zip(url_list, filename_list):

        if download_file(url, filename):
            print(f"File downloaded successfully as {filename}")
        else:
            print("File download failed.")
        
        unzip_gz(filename, filename.replace(".tar", ""))
    

if __name__ == "__main__":
    main()
