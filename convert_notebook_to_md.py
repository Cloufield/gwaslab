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
    
    # Handle text output
    if 'text/plain' in output.get('data', {}):
        text = ''.join(output['data']['text/plain'])
    elif 'text' in output:
        text = ''.join(output['text']) if isinstance(output['text'], list) else str(output['text'])
    elif 'text/html' in output.get('data', {}):
        # For HTML output, we'll just show it as code
        html = ''.join(output['data']['text/html'])
        result.append('```html')
        result.append(html)
        result.append('```')
        return '\n'.join(result)
    else:
        return ''
    
    if not text.strip():
        return ''
    
    # Try to format as table if it looks like one
    formatted = format_table(text)
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
                        md_lines.append(output_md)
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


def get_notebooks_from_mkdocs(mkdocs_path: Path) -> List[Tuple[str, str]]:
    """Extract notebook names from mkdocs.yml Examples section.
    
    Returns list of (title, filename) tuples.
    """
    with open(mkdocs_path, 'r', encoding='utf-8') as f:
        mkdocs = yaml.safe_load(f)
    
    notebooks = []
    nav = mkdocs.get('nav', [])
    
    for section in nav:
        if isinstance(section, dict):
            for key, value in section.items():
                if key == 'Examples' or (isinstance(value, list) and any(
                    isinstance(item, dict) and 'Examples' in item for item in value
                )):
                    # Find Examples section
                    examples_section = None
                    if key == 'Examples':
                        examples_section = value
                    else:
                        for item in value:
                            if isinstance(item, dict) and 'Examples' in item:
                                examples_section = item['Examples']
                                break
                    
                    if examples_section:
                        for item in examples_section:
                            if isinstance(item, dict):
                                for title, filename in item.items():
                                    # Remove .md extension to get base name
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


def convert_from_mkdocs(mkdocs_path: Path, examples_dir: Path, docs_dir: Path):
    """Convert notebooks listed in mkdocs.yml Examples section from examples/ to docs/."""
    notebooks = get_notebooks_from_mkdocs(mkdocs_path)
    
    if not notebooks:
        print("No notebooks found in Examples section of mkdocs.yml")
        return
    
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
  # Convert notebooks from mkdocs.yml Examples section
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
                       help='Convert notebooks listed in mkdocs.yml Examples section from examples/ to docs/')
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
