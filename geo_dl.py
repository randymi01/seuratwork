import tarfile
from concurrent.futures import ThreadPoolExecutor
import os
import click
import requests as re
from tqdm import tqdm
import gzip
import shutil

### TODO
# add single cell data format confirmation and error logging
# output total size


@click.command()
@click.option('-c','--codes', required=True, help='relative file path containing GEO accession codes')
@click.option('-o','--output_path', required=True, help='relative folder path for where data should be stored')
@click.option('-d', '--delete', required=False, default=False, help='delete tar files after extraction')
@click.option('-l', '--download', required=False, default=True, help='Set true to download tar archives')
@click.option('-e', '--extract', required=False, default=False, help='Set true to extract inner tar archives')
def main(codes,output_path, delete, download, extract):
    os.makedirs(output_path, exist_ok=True)

    print("Downloading files...")

    if download:
        with open(codes, "r") as f:
            accession_codes = f.read().splitlines() 
            urls = [geo_url(code) for code in accession_codes]
            tarball = [code+".tar" for code in accession_codes]
            url_files = zip(urls, tarball)
        download_files(url_files, output_path)

    print('\n')
    print("Extracting files...")

    tars = [(output_path+'/'+i, output_path + '/' + i.split('.')[0]) for i in os.listdir(output_path)]


    extract_files(tars, extract)

    print('\n')
    print("Organizing files...")

    folder_paths = [output_path+'/'+i for i in os.listdir(output_path) if os.path.isdir(output_path+'/'+i)]
    
    organize_folders(folder_paths)

    print('\n')
    if delete:
        print("Deleting tar files...")
        for tar, _ in tars:
            os.remove(tar)




def geo_url(code):
    return "https://www.ncbi.nlm.nih.gov/geo/download/?acc="+code+"&format=file"

def download_file(url, output_path):
    try:
        resp = re.get(url, stream = True)
    except:
        print("Error: ", url)

    with open(output_path, "wb") as file, tqdm(
        desc=output_path,
        total=int(resp.headers.get('content-length')),
        unit='iB',
        unit_scale=True,
        unit_divisor=1024,
    ) as pbar:
        for data in resp.iter_content(chunk_size=1024):
            bytes_written = file.write(data)
            pbar.update(bytes_written)

def download_files(url_list, output_path):
    with ThreadPoolExecutor(max_workers=5) as executor:
        for url, ball_name in url_list:
            executor.submit(download_file, url, output_path+'/'+ball_name)


# Extract Tar Files:

def extract_nested_tarfile(outer_tar_path, output_folder, extract_inner_tarfile):
    # Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Open the outer tarfile
    with tarfile.open(outer_tar_path, 'r:*') as outer_tar:
        # Iterate through the members of the outer tarfile
        for member in tqdm(outer_tar.getmembers(), desc = outer_tar_path):
            # Extract the current member
            outer_tar.extract(member, output_folder)
            
            if extract_inner_tarfile:
                inner_tar_path = os.path.join(output_folder, member.name)
                
                # Check if the extracted member is a tarfile
                
                if inner_tar_path.endswith('.gz'):
                    with gzip.open(inner_tar_path) as f_in:
                        with open(inner_tar_path[:-3], 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                else:
                    raise ValueError(f"Unexpected file type: {inner_tar_path}, expected .gz file")
                
                os.remove(inner_tar_path)
            
    # Remove outer tar path            
    # os.remove(outer_tar_path)

# gets the part of the file before barcodes, features (genes), or matrix . filetype
def extract_sample_id(file_name):
    strip_file_name = file_name.split('.')[0]

    if strip_file_name.endswith("barcodes"):
        return strip_file_name[:-9]
    
    if strip_file_name.endswith("features"):
        return strip_file_name[:-9]
    
    if strip_file_name.endswith("genes"):
        return strip_file_name[:-6]
    
    if strip_file_name.endswith("matrix"):
        return strip_file_name[:-7]
    
    return strip_file_name

def extract_seurat_file_name(file_name):
    if file_name.endswith("barcodes.tsv.gz"):
        return "barcodes.tsv.gz"
    
    if file_name.endswith("barcodes.tsv"):
        return "barcodes.tsv"
    
    if file_name.endswith("features.tsv.gz"):
        return "features.tsv.gz"
    
    if file_name.endswith("features.tsv"):
        return "features.tsv"
    
    if file_name.endswith("matrix.mtx.gz"):
        return "matrix.mtx.gz"
    
    if file_name.endswith("matrix.mtx"):
        return "matrix.mtx"
    
    if file_name.endswith("genes.tsv.gz"):
        return "features.tsv.gz"
    
    if file_name.endswith("genes.tsv"):
        return "features.tsv"
    
    return file_name

def organize_folder(folder_path):
    file_struct = {}
    for file in os.listdir(folder_path):
        sample_id = extract_sample_id(file)
        if sample_id not in file_struct:
            file_struct[sample_id] = [file]
        else:
            file_struct[sample_id].append(file)
    # create folder
    for sample_id, files in file_struct.items():
        os.makedirs(folder_path + '/' + sample_id, exist_ok=True)
        for file in files:
            os.rename(folder_path + '/' + file, folder_path + '/' + sample_id + '/' + extract_seurat_file_name(file))

def extract_files(tars, inner_file_extract = False):
    with ThreadPoolExecutor(max_workers=5) as executor:
        for tar_path, output_path in tars:
            executor.submit(extract_nested_tarfile, tar_path, output_path, inner_file_extract)

def organize_folders(folder_paths):
    with ThreadPoolExecutor(max_workers=5) as executor:
        for folder_path in folder_paths:
            executor.submit(organize_folder, folder_path)

if __name__ == "__main__":
    main()





    