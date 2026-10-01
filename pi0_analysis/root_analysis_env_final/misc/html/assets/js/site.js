(() => {
  const body = document.body;
  const theme = localStorage.getItem('nps-theme');
  if (theme === 'light') body.classList.add('light');

  const themeButton = document.querySelector('[data-theme-toggle]');
  if (themeButton) {
    themeButton.textContent = body.classList.contains('light') ? '◐ Dark' : '◑ Light';
    themeButton.addEventListener('click', () => {
      body.classList.toggle('light');
      localStorage.setItem('nps-theme', body.classList.contains('light') ? 'light' : 'dark');
      themeButton.textContent = body.classList.contains('light') ? '◐ Dark' : '◑ Light';
    });
  }

  const menuButton = document.querySelector('[data-menu-toggle]');
  const nav = document.querySelector('.site-nav');
  if (menuButton && nav) menuButton.addEventListener('click', () => nav.classList.toggle('open'));

  const page = body.dataset.page;
  document.querySelectorAll('.site-nav a').forEach(a => {
    if (a.dataset.page === page) a.classList.add('active');
  });

  document.querySelectorAll('pre').forEach(pre => {
    const button = document.createElement('button');
    button.className = 'copy-button'; button.type = 'button'; button.textContent = 'copy';
    button.addEventListener('click', async () => {
      const text = pre.innerText.replace(/^copy\n/, '');
      try { await navigator.clipboard.writeText(text); button.textContent = 'copied'; }
      catch { button.textContent = 'select text'; }
      setTimeout(() => button.textContent = 'copy', 1400);
    });
    pre.appendChild(button);
  });

  const box = document.querySelector('.lightbox');
  const boxImg = box && box.querySelector('img');
  const close = () => box && box.classList.remove('open');
  document.querySelectorAll('[data-lightbox]').forEach(button => button.addEventListener('click', () => {
    if (!box || !boxImg) return;
    const img = button.querySelector('img'); boxImg.src = img.src; boxImg.alt = img.alt; box.classList.add('open');
  }));
  if (box) { box.addEventListener('click', e => { if (e.target === box || e.target.tagName === 'BUTTON') close(); }); }
  document.addEventListener('keydown', e => { if (e.key === 'Escape') close(); });

  document.querySelectorAll('[data-filter]').forEach(button => button.addEventListener('click', () => {
    const filter = button.dataset.filter;
    document.querySelectorAll('[data-filter]').forEach(b => b.classList.toggle('active', b === button));
    document.querySelectorAll('.plot-card[data-category]').forEach(card => {
      card.hidden = filter !== 'all' && card.dataset.category !== filter;
    });
  }));
})();
