import pandas as pd
import requests
from bs4 import BeautifulSoup
import os


def download_file(url, folder):
    response = requests.get(url)
    response.raise_for_status()

    # 从 URL 中提取文件名
    filename = os.path.join(folder, url.split('/')[-1])
    with open(filename, 'wb') as f:
        f.write(response.content)
    print(f"Downloaded: {filename}")


def extract_appendix_data(url, download_folder):
    try:
        response = requests.get(url)
        response.raise_for_status()

        soup = BeautifulSoup(response.content, 'html.parser')

        # 查找附录链接，假设链接在特定标签中
        appendix_links = soup.find_all('a', href=True)

        for link in appendix_links:
            if 'appendix' in link.text.lower():  # 根据文本判断是否为附录文件
                full_link = link['href']
                if not full_link.startswith('http'):
                    full_link = url.rsplit('/', 1)[0] + '/' + full_link  # 处理相对链接
                download_file(full_link, download_folder)

    except Exception as e:
        print(f"Error while fetching data from {url}: {e}")


def main(excel_file, download_folder):
    df = pd.read_excel(excel_file)

    if 'DOI Link' not in df.columns:
        print("The specified column 'DOI Link' does not exist in the Excel file.")
        return

    os.makedirs(download_folder, exist_ok=True)  # 创建下载文件夹

    for doi_url in df['DOI Link']:
        print(f"Processing DOI URL: {doi_url}")
        extract_appendix_data(doi_url, download_folder)


if __name__ == "__main__":
    # 替换为你的 Excel 文件路径和下载文件夹路径
    main("C:\\Users\\andusk\\Desktop\\黄河数据\\test.xlsx", "./downloads")


