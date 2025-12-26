#!/usr/bin/env python3
"""
Convert Jupyter notebooks to Markdown files with proper handling of:
- Code cells with syntax highlighting
- Output text and tables (with preview like pandas)
- Images (saved to docs/images/notebooks/)
- Markdown cells

Usage:
    python convert_notebook_to_md.py --from-mkdocs  # Convert notebooks listed in mkdocs.yml Examples section
    python convert_notebook_to_md.py <notebook.ipynb> [output.md]
    python convert_notebook_to_md.py docs/*.ipynb  # Convert all notebooks in docs/
"""

import json
import sys
import os
import base64
import re
import yaml
import html as html_module
from html.parser import HTMLParser
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import argparse


def extract_image_data(output: Dict[str, Any]) -> Optional[tuple]:
    """Extract image data from notebook output."""
    if 'data' in output:
        data = output['data']
        # Check for PNG
        if 'image/png' in data:
            return ('png', data['image/png'])
        # Check for JPEG
        elif 'image/jpeg' in data:
            return ('jpeg', data['image/jpeg'])
        # Check for SVG
        elif 'image/svg+xml' in data:
            return ('svg', data['image/svg+xml'])
    return None


class TableHTMLParser(HTMLParser):
    """Parser to extract table data from HTML."""
    def __init__(self):
        super().__init__()
        self.headers = []
        self.rows = []
        self.current_row = []
        self.in_table = False
        self.in_thead = False
        self.in_tbody = False
        self.in_tr = False
        self.in_cell = False
        self.current_cell = []
    
    def handle_starttag(self, tag, attrs):
        tag_lower = tag.lower()
        if tag_lower == 'table':
            self.in_table = True
        elif tag_lower == 'thead':
            self.in_thead = True
        elif tag_lower == 'tbody':
            self.in_tbody = True
        elif tag_lower == 'tr':
            self.in_tr = True
            self.current_row = []
        elif tag_lower in ('td', 'th'):
            self.in_cell = True
            self.current_cell = []
    
    def handle_endtag(self, tag):
        tag_lower = tag.lower()
        if tag_lower == 'table':
            self.in_table = False
        elif tag_lower == 'thead':
            self.in_thead = False
        elif tag_lower == 'tbody':
            self.in_tbody = False
        elif tag_lower == 'tr':
            if self.current_row:
                if self.in_thead:
                    # This is a header row
                    if not self.headers:
                        self.headers = self.current_row[:]
                else:
                    # This is a data row
                    self.rows.append(self.current_row[:])
            self.current_row = []
            self.in_tr = False
        elif tag_lower in ('td', 'th'):
            if self.in_cell:
                cell_text = ''.join(self.current_cell).strip()
                self.current_row.append(cell_text)
                self.current_cell = []
                self.in_cell = False
    
    def handle_data(self, data):
        if self.in_cell:
            self.current_cell.append(data)


def parse_html_table(html: str) -> Optional[Tuple[List[str], List[List[str]]]]:
    """Parse HTML table and return (headers, rows)."""
    parser = TableHTMLParser()
    try:
        parser.feed(html)
        if parser.headers and parser.rows:
            return (parser.headers, parser.rows)
    except Exception:
        pass
    return None


def html_table_to_markdown(html: str, max_rows: int = 10, text_plain: Optional[str] = None) -> Optional[str]:
    """Convert HTML table to markdown table format.
    
    Parameters
    ----------
    html : str
        HTML content containing the table
    max_rows : int
        Maximum number of rows to show in preview
    text_plain : str, optional
        Text/plain output that might contain table size info
    """
    # Try to extract table size info from HTML or text/plain (Jupyter format: [N rows x M columns])
    # This should be the original dataframe size, not the number of rows displayed
    table_size_info = None
    # Try multiple patterns to find the size info
    size_patterns = [
        r'\[(\d+)\s+rows?\s+x\s+(\d+)\s+columns?\]',  # Standard format
        r'(\d+)\s+rows?\s+x\s+(\d+)\s+columns?',  # Without brackets
        r'\[(\d+)\s+rows?,\s+(\d+)\s+columns?\]',  # With comma
    ]
    
    # First try HTML
    for pattern in size_patterns:
        match = re.search(pattern, html, re.IGNORECASE)
        if match:
            num_rows = int(match.group(1))
            num_cols = int(match.group(2))
            table_size_info = f"[{num_rows} rows x {num_cols} columns]"
            break
    
    # If not found in HTML, try text/plain output
    if not table_size_info and text_plain:
        for pattern in size_patterns:
            match = re.search(pattern, text_plain, re.IGNORECASE)
            if match:
                num_rows = int(match.group(1))
                num_cols = int(match.group(2))
                table_size_info = f"[{num_rows} rows x {num_cols} columns]"
                break
    
    table_data = parse_html_table(html)
    if not table_data:
        return None
    
    headers, rows = table_data
    
    if not headers or not rows:
        return None
    
    # Detect and remove index column (first column if it looks like an index)
    # Also remove any columns with empty headers (they break markdown parsing)
    remove_first = False
    
    if len(headers) > 1 and len(rows) > 0:
        # Always remove first column if header is empty (definitely an index, breaks parsing)
        first_header = headers[0].strip()
        if not first_header:
            remove_first = True
        # Check if first column contains mostly numeric values (likely index)
        else:
            try:
                # Check more rows to get better sample
                sample_rows = min(50, len(rows))
                first_col_values = [row[0].strip() for row in rows[:sample_rows] if len(row) > 0 and row[0].strip()]
                
                # Also check last few rows if table is large
                if len(rows) > sample_rows:
                    last_col_values = [row[0].strip() for row in rows[-sample_rows:] if len(row) > 0 and row[0].strip()]
                    first_col_values.extend(last_col_values)
                
                if len(first_col_values) >= 5:  # Need at least 5 values to be confident
                    # Check if most values are numeric (integers) - at least 80%
                    numeric_count = sum(1 for val in first_col_values if val.isdigit())
                    if numeric_count >= len(first_col_values) * 0.8:
                        # If mostly numeric, it's likely an index column
                        remove_first = True
            except (ValueError, IndexError):
                pass
    
    # Remove first column if it's an index
    if remove_first and len(headers) > 1:
        headers = headers[1:]
        rows = [row[1:] if len(row) > 1 else row for row in rows]
    
    # Remove any remaining columns with empty headers (from right to left to preserve indices)
    # Empty headers break markdown table parsing in mkdocs
    empty_header_indices = [i for i, header in enumerate(headers) if not header.strip()]
    if empty_header_indices:
        for idx in sorted(empty_header_indices, reverse=True):
            if idx < len(headers):
                headers.pop(idx)
                rows = [row[:idx] + row[idx+1:] if len(row) > idx else row for row in rows]
    
    # Ensure we have at least one header and all headers are non-empty
    if not headers or any(not h.strip() for h in headers):
        # If any header is still empty, filter them out
        valid_indices = [i for i, h in enumerate(headers) if h.strip()]
        if not valid_indices:
            return None
        headers = [headers[i] for i in valid_indices]
        rows = [[row[i] if i < len(row) else '' for i in valid_indices] for row in rows]
    
    # Create markdown table
    md_lines = []
    md_lines.append('| ' + ' | '.join(headers) + ' |')
    md_lines.append('| ' + ' | '.join(['---'] * len(headers)) + ' |')
    
    # Show preview if table is large
    if len(rows) > max_rows:
        # Show first half
        for row in rows[:max_rows//2]:
            # Pad row to match header count
            while len(row) < len(headers):
                row.append('')
            row = row[:len(headers)]
            md_lines.append('| ' + ' | '.join(row) + ' |')
        
        # Add ellipsis row
        ellipsis_row = ['...'] * len(headers)
        md_lines.append('| ' + ' | '.join(ellipsis_row) + ' |')
        
        # Show last half
        for row in rows[-max_rows//2:]:
            # Pad row to match header count
            while len(row) < len(headers):
                row.append('')
            row = row[:len(headers)]
            md_lines.append('| ' + ' | '.join(row) + ' |')
    else:
        # Show all rows
        for row in rows:
            # Pad row to match header count
            while len(row) < len(headers):
                row.append('')
            row = row[:len(headers)]
            md_lines.append('| ' + ' | '.join(row) + ' |')
    
    # Add table size footer (like Jupyter notebooks)
    # Only show if we found the original dataframe size from HTML
    # Don't use parsed row count as it might be a preview
    if table_size_info:
        md_lines.append('')
        md_lines.append(f'*{table_size_info}*')
    
    return '\n'.join(md_lines)


def save_image(image_data: str, image_format: str, notebook_name: str, image_index: int, 
               output_dir: Path, images_dir: Path) -> str:
    """Save image to disk and return markdown image reference."""
    # Create images directory if it doesn't exist
    images_dir.mkdir(parents=True, exist_ok=True)
    
    # Generate image filename
    image_filename = f"{notebook_name}_img_{image_index}.{image_format}"
    image_path = images_dir / image_filename
    
    # Save image
    if image_format == 'svg':
        # SVG is text, not base64
        with open(image_path, 'w', encoding='utf-8') as f:
            f.write(image_data)
    else:
        # PNG/JPEG are base64 encoded
        image_bytes = base64.b64decode(image_data)
        with open(image_path, 'wb') as f:
            f.write(image_bytes)
    
    # Return relative path from output_dir to image
    rel_path = os.path.relpath(image_path, output_dir)
    return rel_path.replace('\\', '/')  # Use forward slashes for markdown


def format_table_preview(text: str, max_rows: int = 10) -> str:
    """Format text output as a markdown table with preview (like pandas).
    
    Shows first max_rows//2 rows, ..., and last max_rows//2 rows if table is large.
    """
    lines = [line for line in text.strip().split('\n') if line.strip()]
    if len(lines) < 2:
        return text
    
    # Check if it looks like a pandas DataFrame output
    # Look for patterns like "[N rows x M columns]" at the end
    pandas_footer_pattern = r'\[(\d+)\s+rows?\s+x\s+\d+\s+columns?\]'
    pandas_footer = None
    data_lines = lines
    
    # Check last few lines for pandas footer
    for i in range(min(3, len(lines))):
        match = re.search(pandas_footer_pattern, lines[-(i+1)], re.IGNORECASE)
        if match:
            pandas_footer = match.group(0)
            num_rows = int(match.group(1))
            # Remove footer from data lines
            data_lines = lines[:-1] if i == 0 else lines[:-(i+1)]
            break
    
    # Check if it looks like a table (has multiple columns separated by spaces/tabs)
    if len(data_lines) < 2:
        return text
    
    first_line = data_lines[0]
    if '\t' in first_line or '  ' in first_line:
        # Try to format as markdown table
        # Split by tabs or multiple spaces
        if '\t' in first_line:
            headers = [h.strip() for h in first_line.split('\t')]
        else:
            # Split by 2+ spaces
            headers = [h.strip() for h in re.split(r'\s{2,}', first_line)]
        
        if len(headers) > 1:
            # Parse data rows
            data_rows = []
            for line in data_lines[1:]:
                if '\t' in line:
                    cells = [c.strip() for c in line.split('\t')]
                else:
                    cells = [c.strip() for c in re.split(r'\s{2,}', line)]
                
                # Pad cells to match header count
                while len(cells) < len(headers):
                    cells.append('')
                cells = cells[:len(headers)]
                
                # Skip empty rows
                if any(c.strip() for c in cells):
                    data_rows.append(cells)
            
            # Create markdown table with preview
            md_lines = []
            md_lines.append('| ' + ' | '.join(headers) + ' |')
            md_lines.append('| ' + ' | '.join(['---'] * len(headers)) + ' |')
            
            # Show preview if table is large
            if len(data_rows) > max_rows and pandas_footer:
                # Show first half
                for row in data_rows[:max_rows//2]:
                    md_lines.append('| ' + ' | '.join(row) + ' |')
                
                # Add ellipsis row
                ellipsis_row = ['...'] * len(headers)
                md_lines.append('| ' + ' | '.join(ellipsis_row) + ' |')
                
                # Show last half
                for row in data_rows[-max_rows//2:]:
                    md_lines.append('| ' + ' | '.join(row) + ' |')
                
                # Add footer if it was a pandas DataFrame
                if pandas_footer:
                    md_lines.append('')
                    md_lines.append(f'*{pandas_footer}*')
            else:
                # Show all rows
                for row in data_rows:
                    md_lines.append('| ' + ' | '.join(row) + ' |')
                
                # Add footer if it was a pandas DataFrame
                if pandas_footer:
                    md_lines.append('')
                    md_lines.append(f'*{pandas_footer}*')
            
            return '\n'.join(md_lines)
    
    return text


def format_table(text: str) -> str:
    """Format text output as a markdown table if it looks like a table (with preview)."""
    return format_table_preview(text, max_rows=10)


def convert_output(output: Dict[str, Any], notebook_name: str, image_index: int,
                   output_dir: Path, images_dir: Path) -> str:
    """Convert notebook output to markdown."""
    result = []
    
    # Check for images
    image_info = extract_image_data(output)
    if image_info:
        image_format, image_data = image_info
        image_path = save_image(image_data, image_format, notebook_name, image_index,
                               output_dir, images_dir)
        result.append(f'![Output image]({image_path})')
        return '\n'.join(result)
    
    # Prioritize HTML output (Jupyter often outputs HTML tables directly)
    if 'text/html' in output.get('data', {}):
        html = ''.join(output['data']['text/html'])
        # Check if it's an HTML table
        if '<table' in html.lower():
            # Get text/plain output if available (for table size info)
            text_plain = None
            if 'text/plain' in output.get('data', {}):
                text_plain = ''.join(output['data']['text/plain'])
            # Parse HTML table and convert to markdown
            markdown_table = html_table_to_markdown(html, max_rows=10, text_plain=text_plain)
            if markdown_table:
                result.append(markdown_table)
                return '\n'.join(result)
            else:
                # If parsing failed, fall back to using HTML directly
                styled_html = f'<div style="font-size: 0.85em; line-height: 1.3;">\n{html}\n</div>'
                result.append(styled_html)
                return '\n'.join(result)
        else:
            # Other HTML, show as code
            result.append('```html')
            result.append(html)
            result.append('```')
            return '\n'.join(result)
    
    # Handle text output
    if 'text/plain' in output.get('data', {}):
        text = ''.join(output['data']['text/plain'])
    elif 'text' in output:
        text = ''.join(output['text']) if isinstance(output['text'], list) else str(output['text'])
    else:
        return ''
    
    if not text.strip():
        return ''
    
    # Try to format as table if it looks like one
    formatted = format_table(text)
    
    # Check if formatted output is a markdown table (contains | characters)
    # If it's a table, don't wrap in code block; otherwise wrap in code block
    is_table = formatted != text and '|' in formatted and formatted.strip().startswith('|')
    
    if is_table:
        # It's a markdown table, add directly without code block
        result.append(formatted)
    else:
        # Regular text output, wrap in code block
        result.append('```')
        result.append(formatted)
        result.append('```')
    
    return '\n'.join(result)


def convert_notebook_to_markdown(notebook_path: Path, output_path: Optional[Path] = None,
                                 images_dir: Optional[Path] = None) -> Path:
    """Convert a Jupyter notebook to a Markdown file."""
    # Read notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        notebook = json.load(f)
    
    # Determine output path
    if output_path is None:
        output_path = notebook_path.with_suffix('.md')
    
    # Determine images directory
    if images_dir is None:
        # Default to docs/images/notebooks/
        docs_dir = output_path.parent
        images_dir = docs_dir / 'images' / 'notebooks'
    
    # Get notebook name for image naming
    notebook_name = notebook_path.stem
    
    # Build markdown content
    md_lines = []
    image_index = 0
    
    for cell in notebook.get('cells', []):
        cell_type = cell.get('cell_type', '')
        source = cell.get('source', [])
        
        if isinstance(source, list):
            cell_content = ''.join(source)
        else:
            cell_content = str(source)
        
        if cell_type == 'markdown':
            # Markdown cell - add directly
            md_lines.append(cell_content)
            md_lines.append('')  # Add blank line
        
        elif cell_type == 'code':
            # Code cell
            if cell_content.strip():
                # Determine language from first line or default to python
                language = 'python'
                if cell_content.strip().startswith('%%'):
                    # Magic commands
                    first_line = cell_content.split('\n')[0]
                    if 'bash' in first_line or 'sh' in first_line:
                        language = 'bash'
                    elif 'r' in first_line or 'R' in first_line:
                        language = 'r'
                
                # Wrap Python code blocks with !!! example for tutorial pages
                if language == 'python':
                    md_lines.append('!!! example')
                    md_lines.append('    ```python')
                    # Indent code content by 4 spaces
                    indented_content = '\n'.join('    ' + line if line else '' for line in cell_content.rstrip().split('\n'))
                    md_lines.append(indented_content)
                    md_lines.append('    ```')
                else:
                    md_lines.append(f'```{language}')
                    md_lines.append(cell_content.rstrip())
                    md_lines.append('```')
                md_lines.append('')
            
            # Handle outputs
            outputs = cell.get('outputs', [])
            for output in outputs:
                output_type = output.get('output_type', '')
                
                if output_type == 'stream':
                    # Stream output (stdout/stderr)
                    name = output.get('name', 'stdout')
                    text = ''.join(output.get('text', []))
                    if text.strip():
                        md_lines.append(f'**{name}:**')
                        md_lines.append('```')
                        md_lines.append(text.rstrip())
                        md_lines.append('```')
                        md_lines.append('')
                
                elif output_type == 'execute_result' or output_type == 'display_data':
                    # Execute result or display data
                    output_md = convert_output(output, notebook_name, image_index,
                                              output_path.parent, images_dir)
                    if output_md:
                        # Split by newlines to ensure each line is separate
                        # This is important for markdown tables to render correctly
                        output_lines = output_md.split('\n')
                        md_lines.extend(output_lines)
                        md_lines.append('')
                        if 'image' in output_md:
                            image_index += 1
    
    # Write markdown file
    md_content = '\n'.join(md_lines)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(md_content)
    
    print(f"Converted {notebook_path} -> {output_path}")
    if image_index > 0:
        print(f"  Saved {image_index} image(s) to {images_dir}")
    
    return output_path


def get_notebooks_from_mkdocs(mkdocs_path: Path, sections: List[str] = None) -> List[Tuple[str, str]]:
    """Extract notebook names from mkdocs.yml sections.
    
    Parameters
    ----------
    mkdocs_path : Path
        Path to mkdocs.yml file
    sections : List[str], optional
        List of section names to extract from. If None, extracts from all sections.
        Default: ['Examples', 'Start'] to include Examples and Start sections.
    
    Returns
    -------
    List[Tuple[str, str]]
        List of (title, filename) tuples.
    """
    if sections is None:
        sections = ['Examples', 'Start']
    
    with open(mkdocs_path, 'r', encoding='utf-8') as f:
        mkdocs = yaml.safe_load(f)
    
    notebooks = []
    nav = mkdocs.get('nav', [])
    
    for section in nav:
        if isinstance(section, dict):
            for key, value in section.items():
                # Check if this section should be processed
                if key in sections:
                    # Process items in this section
                    if isinstance(value, list):
                        for item in value:
                            if isinstance(item, dict):
                                for title, filename in item.items():
                                    # Remove .md extension to get base name
                                    base_name = filename.replace('.md', '')
                                    notebooks.append((title, base_name))
                    elif isinstance(value, str):
                        # Direct string value (e.g., "Tutorial: tutorial_v4.md")
                        # This format is less common but handle it
                        if value.endswith('.md'):
                            base_name = value.replace('.md', '')
                            notebooks.append((key, base_name))
                elif isinstance(value, list) and any(
                    isinstance(item, dict) and section_name in item 
                    for item in value for section_name in sections
                ):
                    # Handle nested sections (e.g., when Examples is nested)
                    for item in value:
                        if isinstance(item, dict):
                            for section_name in sections:
                                if section_name in item:
                                    section_items = item[section_name]
                                    if isinstance(section_items, list):
                                        for sub_item in section_items:
                                            if isinstance(sub_item, dict):
                                                for title, filename in sub_item.items():
                                                    base_name = filename.replace('.md', '')
                                                    notebooks.append((title, base_name))
    
    return notebooks


def find_notebook_in_examples(notebook_name: str, examples_dir: Path) -> Optional[Path]:
    """Find a notebook in the examples directory."""
    # Try exact match first
    for ipynb_file in examples_dir.rglob(f"{notebook_name}.ipynb"):
        # Skip checkpoint files
        if '.ipynb_checkpoints' not in str(ipynb_file):
            return ipynb_file
    
    return None


def convert_from_mkdocs(mkdocs_path: Path, examples_dir: Path, docs_dir: Path, sections: List[str] = None):
    """Convert notebooks listed in mkdocs.yml sections from examples/ to docs/.
    
    Parameters
    ----------
    mkdocs_path : Path
        Path to mkdocs.yml file
    examples_dir : Path
        Directory containing notebook files
    docs_dir : Path
        Directory to output markdown files
    sections : List[str], optional
        List of section names to process. If None, processes 'Examples' and 'Start' sections.
    """
    notebooks = get_notebooks_from_mkdocs(mkdocs_path, sections=sections)
    
    if not notebooks:
        section_names = ', '.join(sections) if sections else 'specified'
        print(f"No notebooks found in {section_names} section(s) of mkdocs.yml")
        return
    
    section_names = ', '.join(sections) if sections else 'specified'
    print(f"Found {len(notebooks)} notebooks in {section_names} section(s)")
    
    print(f"Found {len(notebooks)} notebooks in Examples section")
    
    images_dir = docs_dir / 'images' / 'notebooks'
    
    for title, notebook_name in notebooks:
        print(f"\nProcessing: {title} ({notebook_name})")
        
        # Find notebook in examples directory
        notebook_path = find_notebook_in_examples(notebook_name, examples_dir)
        
        if not notebook_path:
            print(f"  Warning: Could not find {notebook_name}.ipynb in {examples_dir}")
            continue
        
        # Output to docs directory
        output_path = docs_dir / f"{notebook_name}.md"
        
        convert_notebook_to_markdown(notebook_path, output_path, images_dir)


def main():
    parser = argparse.ArgumentParser(
        description='Convert Jupyter notebooks to Markdown files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert notebooks from mkdocs.yml Examples and Start sections
  python convert_notebook_to_md.py --from-mkdocs
  
  # Convert a single notebook
  python convert_notebook_to_md.py docs/standardization_workflow.ipynb
  
  # Convert with custom output path
  python convert_notebook_to_md.py notebook.ipynb output.md
  
  # Convert all notebooks in docs/
  python convert_notebook_to_md.py docs/*.ipynb
        """
    )
    parser.add_argument('notebooks', nargs='*', help='Notebook file(s) to convert')
    parser.add_argument('-o', '--output', help='Output markdown file (only for single notebook)')
    parser.add_argument('--images-dir', help='Directory to save images (default: docs/images/notebooks/)')
    parser.add_argument('--from-mkdocs', action='store_true',
                       help='Convert notebooks listed in mkdocs.yml Examples and Start sections from examples/ to docs/')
    parser.add_argument('--mkdocs-path', default='mkdocs.yml',
                       help='Path to mkdocs.yml (default: mkdocs.yml)')
    parser.add_argument('--examples-dir', default='examples',
                       help='Path to examples directory (default: examples)')
    parser.add_argument('--docs-dir', default='docs',
                       help='Path to docs directory (default: docs)')
    
    args = parser.parse_args()
    
    # Handle --from-mkdocs mode
    if args.from_mkdocs:
        script_dir = Path(__file__).parent
        mkdocs_path = script_dir / args.mkdocs_path
        examples_dir = script_dir / args.examples_dir
        docs_dir = script_dir / args.docs_dir
        
        if not mkdocs_path.exists():
            print(f"Error: {mkdocs_path} does not exist")
            sys.exit(1)
        
        convert_from_mkdocs(mkdocs_path, examples_dir, docs_dir)
        return
    
    # Handle regular conversion mode
    if not args.notebooks:
        parser.print_help()
        sys.exit(1)
    
    # Handle images directory
    images_dir = None
    if args.images_dir:
        images_dir = Path(args.images_dir)
    
    # Convert notebooks
    if len(args.notebooks) == 1 and args.output:
        # Single notebook with custom output
        notebook_path = Path(args.notebooks[0])
        output_path = Path(args.output)
        convert_notebook_to_markdown(notebook_path, output_path, images_dir)
    else:
        # Multiple notebooks or default output
        for notebook_arg in args.notebooks:
            notebook_path = Path(notebook_arg)
            if not notebook_path.exists():
                print(f"Warning: {notebook_path} does not exist, skipping...")
                continue
            
            output_path = None if args.output and len(args.notebooks) == 1 else None
            convert_notebook_to_markdown(notebook_path, output_path, images_dir)


if __name__ == '__main__':
    main()
