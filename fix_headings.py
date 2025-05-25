import nbformat

input_path = "notebooks/project_report.ipynb"
output_path = "notebooks/project_report.ipynb"  # overwrite the original

notebook = nbformat.read(open(input_path, "r", encoding="utf-8"), as_version=4)

for cell in notebook.cells:
    if cell.cell_type == "markdown":
        lines = cell.source.splitlines()
        new_lines = []
        for line in lines:
            if line.startswith("### "):  # Promote from H3 to H2
                new_lines.append("## " + line[4:])
            else:
                new_lines.append(line)
        cell.source = "\n".join(new_lines)

nbformat.write(notebook, open(output_path, "w", encoding="utf-8"))
print("✅ Headings adjusted for better GitHub navigation.")
