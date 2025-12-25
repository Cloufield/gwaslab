#!/usr/bin/env python3
"""
Convert plain code blocks (```) to Python code blocks (```python) in markdown files.
"""

import re
import sys
from pathlib import Path

def convert_code_blocks(content: str) -> str:
    """
    Convert plain code blocks (```) to Python code blocks (```python).
    Only converts opening delimiters that don't already have a language specified.
    """
    lines = content.split('\n')
    result = []
    in_code_block = False
    i = 0
    
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        
        # Check if this is a code block delimiter
        if stripped.startswith('```'):
            # Check if it has a language specified (e.g., ```python, ```bash)
            if len(stripped) > 3 and stripped[3:].strip() and not stripped[3:].strip().startswith('`'):
                # Already has a language, keep as is
                result.append(line)
                in_code_block = not in_code_block
            elif stripped == '```':
                # Plain code block delimiter
                if in_code_block:
                    # This is a closing delimiter, keep as is
                    result.append(line)
                    in_code_block = False
                else:
                    # This is an opening delimiter, convert to python
                    # Preserve indentation
                    indent = len(line) - len(line.lstrip())
                    result.append(' ' * indent + '```python')
                    in_code_block = True
            else:
                # Multiple backticks or other case, keep as is
                result.append(line)
        else:
            result.append(line)
        
        i += 1
    
    return '\n'.join(result)

def process_file(file_path: Path) -> bool:
    """Process a single markdown file."""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        original_content = content
        new_content = convert_code_blocks(content)
        
        if original_content != new_content:
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(new_content)
            return True
        return False
    except Exception as e:
        print(f"Error processing {file_path}: {e}", file=sys.stderr)
        return False

def main():
    """Main function to process all markdown files in docs/."""
    docs_dir = Path('docs')
    
    if not docs_dir.exists():
        print(f"Error: {docs_dir} directory not found", file=sys.stderr)
        sys.exit(1)
    
    # Find all markdown files
    md_files = list(docs_dir.glob('*.md'))
    
    if not md_files:
        print("No markdown files found in docs/")
        return
    
    print(f"Found {len(md_files)} markdown files")
    
    modified_count = 0
    for md_file in sorted(md_files):
        if process_file(md_file):
            print(f"Modified: {md_file}")
            modified_count += 1
        else:
            print(f"No changes: {md_file}")
    
    print(f"\nTotal files modified: {modified_count}/{len(md_files)}")

if __name__ == '__main__':
    main()

