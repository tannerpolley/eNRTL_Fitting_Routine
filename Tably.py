import csv
import os


PREAMBLE = r"""\documentclass[margin=.1cm]{standalone}
\usepackage[labelformat=empty]{caption}
\usepackage{booktabs}
\captionsetup{font=small}
\begin{document}
\begin{minipage}{\linewidth}
\centering
"""

HEADER = r"""
\captionof{table}{caption}
\scriptsize
{indent}\begin{{tabular}}{{@{{}}{align}@{{}}}}
{indent}{indent}\toprule"""

FOOTER = r"""{indent}{indent}\bottomrule
{indent}\end{{tabular}}
\bigskip
\end{minipage}
"""

LABEL = '\n{indent}\\label{{{label}}}'
CAPTION = '{{{caption}}}'


class Tably:
    """Object which holds parsed arguments.

    Methods:
        run: selects the appropriate methods to generate LaTeX code/files
        create_table: for each specified file, creates a LaTeX table
        create_row: creates a row based on `line` content
        combine_tables: combines all tables from input files together
        save_single_table: creates and saves a single LaTeX table
        get_units: writes the units as a row of the LaTeX table
    """

    def __init__(self, files=['Parameters_table.csv'], title='title', folder=os.getcwd()):

        self.files = files
        self.no_header = False
        self.caption = title
        self.label = None
        self.align = 'clrrr'
        self.no_indent = False
        self.outfile = None
        self.separate_outfiles = folder
        self.skip = 0
        self.preamble = True
        self.sep = get_sep(',')
        self.units = None
        self.fragment = False
        self.fragment_skip_header = False
        self.replace = False
        self.tex_str = escape if not False else lambda x: x

    def run(self):
        """The main method.

        If all tables need to be put into a single file,
        calls `combine_tables` method to generate LaTeX code
        and then calls `save_content` function if `outfile` is provided;
        otherwise, prints to the console.
        If each table needs to be put into a separate file,
        calls `save_single_table` method to create and save each table separately.
        """

        if self.fragment_skip_header:
            self.skip = 1
            self.no_header = True
            self.fragment = True

        if self.fragment:
            self.no_indent = True
            self.label = None
            self.preamble = False

        # if all tables need to be put into one file
        if self.outfile or self.separate_outfiles is None:
            final_content = self.combine_tables()
            if not final_content:
                return
            if self.outfile:
                try:
                    save_content(final_content, self.outfile, self.replace)
                except FileNotFoundError:
                    print('{} is not a valid/known path. Could not save there.'.format(self.outfile))
            else:
                return final_content

        # if -oo is passed (could be [])
        if self.separate_outfiles is not None:
            outs = self.separate_outfiles
            if len(outs) == 0:
                outs = [os.path.splitext(file)[0]+'.tex' for file in self.files ]
            elif os.path.isdir(outs[0]):
                outs = [os.path.join(outs[0], os.path.splitext(os.path.basename(file))[0])+'.tex' for file in self.files]
            elif len(outs) != len(self.files):
                print('WARNING: Number of .csv files and number of output files do not match!')
            for file, out in zip(self.files, outs):
                self.save_single_table(file, out)

    def create_table(self, file):
        """Creates a table from a given .csv file.

        This method gives the procedure of converting a .csv file to a LaTeX table.
        Unless -f is specified, the output is a ready-to-use LaTeX table environment.
        All other methods that need to obtain a LaTeX table from a .csv file call this method.
        """
        rows = []
        indent = 4*' ' if not self.no_indent else ''

        try:
            with open(file) as infile:
                for i, columns in enumerate(csv.reader(infile, delimiter=self.sep)):
                    if i < self.skip:
                        continue
                    rows.append(self.create_row(columns, indent))
        except FileNotFoundError:
            print("File {} doesn't exist!!\n".format(file))
            return ''
        if not rows:
            print("No table created from the {} file. Check if the file is empty "
                  "or you used too high skip value.\n".format(file))
            return ''

        if not self.no_header:
            rows.insert(1, r'{0}{0}\midrule'.format(indent))
            if self.units:
                rows[0] = rows[0] + r'\relax' # fixes problem with \[
                units = self.get_units()
                rows.insert(1, r'{0}{0}{1} \\'.format(indent, units))

        content = '\n'.join(rows)
        if not self.fragment:
            header = HEADER.format(
                table='{table}',
                label=add_label(self.label, indent),
                caption=add_caption(self.caption, indent),
                align=format_alignment(self.align, len(columns)),
                indent=indent,
                )
            footer = FOOTER.format(indent=indent, minipage='{minipage}')
            return '\n'.join((header, content, footer))
        else:
            return content

    def create_row(self, line, indent):
        """Creates a row based on `line` content"""
        return r'{indent}{indent}{content} \\'.format(
             indent=indent,
             content=' & '.join(self.tex_str(line)))

    def combine_tables(self):
        """Combine all tables together and add a preamble if required.

        Unless -oo is specified, this is how input tables are arranged.
        """
        all_tables = []
        if self.label and len(self.files) > 1:
            all_tables.append("% don't forget to manually re-label the tables")

        for file in self.files:
            table = self.create_table(file)
            if table:
                all_tables.append(table)
        if not all_tables:
            return None
        if self.preamble:
            all_tables.insert(0, PREAMBLE)
            all_tables.append('\\end{document}\n')
        return '\n\n'.join(all_tables)

    def save_single_table(self, file, out):
        """Creates and saves a single LaTeX table"""
        table = [self.create_table(file)]
        if table:
            if self.preamble:
                table.insert(0, PREAMBLE)
                table.append('\\end{document}\n')
            final_content = '\n\n'.join(table)
            try:
                save_content(final_content, out, self.replace)
            except FileNotFoundError:
                print('{} is not a valid/known path. Could not save there.'.format(out))

    def get_units(self):
        """Writes the units as a row of the LaTeX table"""
        formatted_units = []
        for unit in self.tex_str(self.units):
            if unit in '-/0':
                formatted_units.append('')
            else:
                formatted_units.append('[{}]'.format(unit))
        return ' & '.join(formatted_units)


def get_sep(sep):
    if sep.lower() in ['t', 'tab', '\\t']:
        return '\t'
    elif sep.lower() in ['s', 'semi', ';']:
        return ';'
    elif sep.lower() in ['c', 'comma', ',']:
        return ','
    else:
        return sep


def escape(line):
    """Escapes special LaTeX characters by prefixing them with backslash"""
    for char in '#%&':
        line = [column.replace(char, '\\'+char) for column in line]
    return line


def format_alignment(align, length):
    """Makes sure that provided alignment is valid:
    1. the length of alignment is either 1 or the same as the number of columns
    2. valid characters are `l`, `c` and `r`

    If there is an invalid character, all columns are set to centered alignment.
    If alignment length is too long, it is stripped to fit the number of columns.
    If alignment length is too short, it is padded with `c` for the missing
    columns.
    """
    if any(ch not in 'lcr' for ch in align):
        align = 'c'

    if len(align) == 1:
        return length * align
    elif len(align) == length:
        return align
    else:
        return '{:c<{l}.{l}}'.format(align, l=length)


def add_label(label, indent):
    """Creates a table label"""
    return LABEL.format(label=label, indent=indent) if label else ''


def add_caption(caption, indent):
    """Creates a table caption"""
    return CAPTION.format(caption=caption, indent=indent) if caption else ''


def save_content(content, outfile, replace):
    """Saves the content to a file.

    If an existing file is provided, the content is appended to the end
    of the file by default. If -r is passed, the file is overwritten.
    """
    if replace:
        with open(outfile, 'w') as out:
            out.writelines(content)
        print('The content is written to', outfile)
    else:
        with open(outfile, 'a') as out:
            out.writelines(content)
        # print('The content is appended to', outfile)


def get_latex(files, title, folder):
    tably = Tably(files=files, title=title, folder=folder)
    return tably.run()

