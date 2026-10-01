// Send the visible input strings with the click, avoiding FloatText comm timing.
export default {render({model, el}) {
  const button = document.createElement('button');
  button.type = 'button';
  button.className = 'jupyter-button widget-button';
  button.textContent = 'Apply outer limits';
  button.style.cssText = 'width:180px;flex:0 0 180px';
  const status = document.createElement('span');
  status.setAttribute('role', 'status');
  status.style.cssText = 'font:12px sans-serif;color:#333;max-width:190px;overflow-wrap:anywhere';
  el.style.cssText = 'display:flex;align-items:center;gap:6px';
  el.append(button, status);
  let pending = false;
  const apply = () => {
    if (pending) return;
    const root = el.closest('.lt-dashboard');
    const limits = {};
    for (const [key, prefix] of [['tprime', 't'], ['Q2', 'q'], ['xB', 'x']]) {
      const fields = ['min', 'max'].map(bound =>
        root?.querySelector(`.outer-limit-${prefix}-${bound} input`));
      if (fields.some(field => !field)) {
        status.textContent = 'Outer limit fields are unavailable.';
        return;
      }
      limits[key] = fields.map(field => field.value);
    }
    pending = true;
    button.disabled = true;
    status.textContent = 'Applying...';
    model.send({kind: 'apply_limits', limits});
  };
  const acknowledge = content => {
    if (content.kind !== 'limits_ack') return;
    pending = false;
    button.disabled = false;
    status.textContent = content.message;
  };
  button.addEventListener('click', apply);
  model.on('msg:custom', acknowledge);
  return () => {
    button.removeEventListener('click', apply);
    model.off('msg:custom', acknowledge);
  };
}};
