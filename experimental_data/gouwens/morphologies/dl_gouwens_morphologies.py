import os
import requests
from urllib.parse import urljoin
from bs4 import BeautifulSoup

def ensure_directory_exists(directory):
    """
    Ensure that the specified directory exists. If not, create it.

    :param directory: Path of the directory.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

def download_directory(url, target_directory):
    """
    Download all files from the specified directory URL and save them in the target directory.

    :param url: URL of the directory to download.
    :param target_directory: Path of the directory to save the downloaded files.
    """
    try:
        # Ensure target directory exists
        ensure_directory_exists(target_directory)

        # Send HTTP GET request to the directory URL
        response = requests.get(url)
        response.raise_for_status()  # Raise an exception for HTTP errors

        # Parse HTML content using BeautifulSoup
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find all anchor tags (links)
        for link in soup.find_all('a'):
            href = link.get('href').strip()
            if href == '../':
                continue  # Skip parent directory link
            file_url = urljoin(url, href)
            file_name = os.path.basename(file_url)  # Extract file name from URL
            file_path = os.path.join(target_directory, file_name)

            # Download each file from the directory
            download_file(file_url, file_path)

    except requests.RequestException as e:
        print(f"Error downloading directory from {url}: {e}")

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

# def download_cell_morphology(cell_class, cell_name):
#     """
#     Download the morphology file for a given cell from the specified cell class.

#     :param cell_class: The class of the cell (e.g., 'excitatory' or 'inhibitory').
#     :param cell_name: The name of the cell.
#     """
#     path = f"{cell_class}/{cell_name}.SWC"
#     document_url = f"https://download.brainimagelibrary.org/3a/88/3a88a7687ab66069/{path}"
#     filename = path
#     download_file(document_url, filename)

def main():
    """
    Main function to download cell morphology files.
    """
    # # Read metadata from a CSV file
    # meta_data = pd.read_csv('../m1_patchseq_meta_data.csv', index_col='Cell')

    # create_directory("excitatory")
    # create_directory("excitatory")

    # for cell_name in meta_data.index:
    #     print(f"Downloading morphology file for cell {cell_name}...")
    #     try:
    #         download_cell_morphology('excitatory', cell_name)
    #         print('Downloaded excitatory morphology')
    #     except Exception as e:
    #         print(f"Failed to download excitatory morphology: {e}")
    #         try:
    #             download_cell_morphology('inhibitory', cell_name)
    #             print('Downloaded inhibitory morphology')
    #         except Exception as e:
    #             print(f"Failed to download inhibitory morphology: {e}")
    #     print('__________')

    directory_url_list = [
        "https://download.brainimagelibrary.org/biccn/zeng/pseq/morph/200526/"
    ]

    target_directory_list = [
        "inhibitory_files"
    ]

    for directory_url, target_directory in zip(directory_url_list, target_directory_list):
        download_directory(directory_url, target_directory)

if __name__ == "__main__":
    main()
