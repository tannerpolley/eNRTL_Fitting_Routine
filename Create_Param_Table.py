import pandas as pd
import os
from pdf2image import convert_from_path
from Tably import get_latex
from pdflatex import PDFLaTeX
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np


def create_param_table(df, title='Table Title', folder=r'data\Parameters', plot=True):

    def insert_file(file_name):
        return folder + r'\\' + file_name
    df_table = df.drop(columns=['Object', 'Object_Name', ])

    df_1 = df_table[(df_table['Name'] != "$k_{car,4}$") & (df_table['Name'] != "$k_{bic,4}$")]
    reaction_4 = [col for col in df_table['Name'] if '4' in col]

    if len(reaction_4) > 0:
        df_2 = df_table[(df_table['Name'].isin(reaction_4))]

    df_table['Value'] = df_table['Value'].astype('str')
    df_table['Uncertainty'] = df_table['Uncertainty'].astype('str')
    df_table.loc[df_1.index, 'Value'] = df_1['Value'].map('{:.1f}'.format)
    df_table.loc[df_1.index, 'Uncertainty'] = df_1['Uncertainty'].map('{:.1f}'.format)

    if len(reaction_4) > 0:
        df_table.loc[df_2.index, 'Value'] = df_2['Value'].map('{:.4f}'.format)
        df_table.loc[df_2.index, 'Uncertainty'] = df_2['Uncertainty'].map('{:.4f}'.format)

    for i, row in df_table.iterrows():
        df_table.loc[i, 'Uncertainty'] = '$\pm$' + ' ' + row["Uncertainty"]

    df_table.sort_values(by=['Name'], ascending=True, inplace=True)
    df_table['Percent'] = df_table['Percent'].map('{:.0%}'.format)

    df_table.to_csv(insert_file('Parameters_table.csv'), index=False)

    # Gets a string of latex code that creates a table from a csv file and stores in to a tex file
    get_latex([insert_file('Parameters_table.csv')], title, [folder])
    os.remove(insert_file('Parameters_table.csv'))

    # Runs the latex file and writes the bytes output to a pdf file
    pdfl = PDFLaTeX.from_texfile(insert_file('Parameters_table.tex'))
    pdfl.set_output_directory(folder)
    os.remove(insert_file('Parameters_table.tex'))

    pdfl.set_interaction_mode()
    pdf, log, completed_process = pdfl.create_pdf(keep_pdf_file=True, keep_log_file=False)
    with open(insert_file('Parameters_table.pdf'), 'wb') as f:
        f.write(pdf)

    # Converts the pdf into an image
    image = convert_from_path(insert_file('Parameters_table.pdf'), dpi=200)
    os.remove(insert_file('Parameters_table.pdf'))

    # Saves the image object as an image file
    image[0].save(insert_file('Parameters.png'))
    print('New Parameter Table saved')
    img = np.asarray(Image.open(insert_file('Parameters.png')))

    if plot:
        plt.figure(figsize=(20, 20))
        plt.imshow(img)
        plt.axis('off')
        plt.tight_layout()
        plt.show()

if __name__ == '__main__':
    df_fit = pd.read_csv(r'data\Parameters\Parameters_fit.csv')
    create_param_table(df_fit, 'Testing title')
