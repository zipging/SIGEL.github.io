document.addEventListener('DOMContentLoaded', function() {
    new ClipboardJS('.btn'); // 确保按钮的类名与下面的 HTML 对应
  
    document.querySelectorAll('pre code').forEach(function(block) {
      var btn = document.createElement('button');
      btn.className = 'btn';
      btn.textContent = 'Copy';
      btn.setAttribute('data-clipboard-target', '#'+block.id);
      block.parentNode.insertBefore(btn, block);
    });
  });
  