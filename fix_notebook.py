import nbformat

input_path = "notebooks/project_report.ipynb"
output_path = "notebooks/project_report_cleaned_markdown.ipynb"

# Load notebook
with open(input_path, "r", encoding="utf-8") as f:
    notebook = nbformat.read(f, as_version=4)

# Convert code cells to markdown
for cell in notebook.cells:
    if cell.cell_type == "code" and cell.source.strip():
        source = cell.source.strip()
        cell.cell_type = "markdown"
        cell.source = f"```python\n{source}\n```"
        cell.outputs = []
        cell.execution_count = None

# Save new notebook
with open(output_path, "w", encoding="utf-8") as f:
    nbformat.write(notebook, f)

print(f"✅ Cleaned: {output_path}")
