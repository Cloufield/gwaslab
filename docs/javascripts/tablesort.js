// Initialize TableSort only on actual markdown tables, not code blocks
(function() {
    'use strict';
    
    // Wait for DOM to be ready and code highlighting to complete
    function initTableSort() {
        // Only target tables that are NOT inside code blocks
        // Select all tables in the document
        const allTables = document.querySelectorAll('table');
        
        // Filter out tables that are inside <pre> or <code> elements (code blocks)
        const markdownTables = Array.from(allTables).filter(table => {
            // Check if table is inside a code block
            // Look for common code block containers
            const isInCodeBlock = table.closest('pre') || 
                                 table.closest('code') ||
                                 table.closest('.highlight') ||
                                 table.closest('.codehilite') ||
                                 table.closest('.highlighttable') ||
                                 table.closest('[class*="highlight"]');
            
            if (isInCodeBlock) {
                return false; // Skip tables inside code blocks
            }
            
            // Additional check: if table is inside a <pre> tag (common for code blocks)
            let parent = table.parentElement;
            while (parent && parent !== document.body) {
                if (parent.tagName === 'PRE' || 
                    parent.tagName === 'CODE' ||
                    (parent.classList && (
                        parent.classList.contains('highlight') ||
                        parent.classList.contains('codehilite') ||
                        Array.from(parent.classList).some(cls => cls.includes('highlight'))
                    ))) {
                    return false;
                }
                parent = parent.parentElement;
            }
            return true; // This is a real markdown table
        });
        
        // Initialize TableSort only on markdown tables
        if (typeof TableSort !== 'undefined' && markdownTables.length > 0) {
            markdownTables.forEach(table => {
                try {
                    // Only initialize if not already initialized
                    if (!table.hasAttribute('data-tablesort-initialized')) {
                        new TableSort(table);
                        table.setAttribute('data-tablesort-initialized', 'true');
                    }
                } catch (e) {
                    console.warn('Failed to initialize TableSort on table:', e);
                }
            });
        }
    }
    
    // Run when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', function() {
            // Wait a bit for code highlighting to complete
            setTimeout(initTableSort, 100);
        });
    } else {
        // DOM already loaded, wait for highlighting
        setTimeout(initTableSort, 100);
    }
    
    // Also run after a longer delay to catch any dynamically loaded content
    // But ensure this happens after code highlighting is complete
    setTimeout(initTableSort, 1000);
})();

