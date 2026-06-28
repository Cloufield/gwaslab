// Protect code blocks from being modified by other scripts
(function() {
    'use strict';
    
    // Function to protect code blocks
    function protectCodeBlocks() {
        // Find all code blocks
        const codeBlocks = document.querySelectorAll('pre code, .highlight code, .codehilite code');
        
        codeBlocks.forEach(code => {
            // Ensure white-space is preserved
            if (code.style.whiteSpace !== 'pre') {
                code.style.whiteSpace = 'pre';
            }
            
            // Ensure word-break is normal
            if (code.style.wordBreak !== 'normal') {
                code.style.wordBreak = 'normal';
            }
            
            // Ensure overflow-wrap is normal
            if (code.style.overflowWrap !== 'normal') {
                code.style.overflowWrap = 'normal';
            }
            
            // Ensure word-wrap is normal
            if (code.style.wordWrap !== 'normal') {
                code.style.wordWrap = 'normal';
            }
            
            // Prevent hyphens
            code.style.hyphens = 'none';
            code.style.webkitHyphens = 'none';
            code.style.mozHyphens = 'none';
            code.style.msHyphens = 'none';
        });
        
        // Also protect pre elements
        const preElements = document.querySelectorAll('pre, .highlight pre, .codehilite pre');
        preElements.forEach(pre => {
            if (pre.style.whiteSpace !== 'pre') {
                pre.style.whiteSpace = 'pre';
            }
            if (pre.style.wordBreak !== 'normal') {
                pre.style.wordBreak = 'normal';
            }
            if (pre.style.overflowWrap !== 'normal') {
                pre.style.overflowWrap = 'normal';
            }
        });
    }
    
    // Run immediately
    protectCodeBlocks();
    
    // Run after DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', protectCodeBlocks);
    } else {
        protectCodeBlocks();
    }
    
    // Run after a delay to catch any scripts that modify code blocks
    setTimeout(protectCodeBlocks, 100);
    setTimeout(protectCodeBlocks, 500);
    setTimeout(protectCodeBlocks, 1000);
    setTimeout(protectCodeBlocks, 2000);
    
    // Watch for changes to code blocks
    const observer = new MutationObserver(function(mutations) {
        let shouldProtect = false;
        mutations.forEach(function(mutation) {
            if (mutation.type === 'attributes' || mutation.type === 'childList') {
                const target = mutation.target;
                if (target.tagName === 'CODE' || target.tagName === 'PRE' || 
                    target.classList.contains('highlight') || 
                    target.classList.contains('codehilite')) {
                    shouldProtect = true;
                }
                // Check if any added nodes are code blocks
                if (mutation.addedNodes) {
                    mutation.addedNodes.forEach(function(node) {
                        if (node.nodeType === 1 && (
                            node.tagName === 'CODE' || 
                            node.tagName === 'PRE' ||
                            node.classList.contains('highlight') ||
                            node.classList.contains('codehilite')
                        )) {
                            shouldProtect = true;
                        }
                    });
                }
            }
        });
        
        if (shouldProtect) {
            setTimeout(protectCodeBlocks, 10);
        }
    });
    
    // Observe the document for changes
    observer.observe(document.body, {
        childList: true,
        subtree: true,
        attributes: true,
        attributeFilter: ['style', 'class']
    });
})();

