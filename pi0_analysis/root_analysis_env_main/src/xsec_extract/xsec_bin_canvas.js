// Exact axis coordinates via the SVG transform, independent of data traces.
export default {render({model, el}) {
  el.classList.add('lt-canvas');
  const message = document.createElement('div');
  message.setAttribute('role', 'status');
  message.style.cssText = 'min-height:22px;font:13px sans-serif;color:#333';
  const host = document.createElement('div');
  host.style.cssText = 'width:100%;max-width:1100px';
  el.append(message, host);
  let pending = false;
  const queue = [];
  const sendNext = () => {
    if (pending || !queue.length) return;
    pending = true;
    host.style.opacity = '0.8';
    model.send(queue.shift());
  };
  const draw = () => {
    host.innerHTML = model.get('svg');
    const svg = host.querySelector('svg');
    if (!svg) return;
    svg.style.cssText = 'width:100%;height:auto;display:block;touch-action:manipulation';
    svg.setAttribute('aria-label', 'Linked LT bin selection plots');
    // No hover elements. Only the two axes accept pointer clicks.
    for (const axis of model.get('geometry')) {
      const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
      for (const [k, v] of Object.entries({x:axis.left, y:axis.top,
          width:axis.width, height:axis.height, fill:'transparent'})) rect.setAttribute(k, v);
      rect.setAttribute('data-view', axis.name);
      rect.style.cursor = 'crosshair';
      rect.addEventListener('click', event => {
        const point = new DOMPoint(event.clientX, event.clientY).matrixTransform(svg.getScreenCTM().inverse());
        const x = axis.xlim[0] + (point.x-axis.left)/axis.width*(axis.xlim[1]-axis.xlim[0]);
        const y = axis.ylim[1] - (point.y-axis.top)/axis.height*(axis.ylim[1]-axis.ylim[0]);
        queue.push({kind:'place', view:axis.name, x, y});
        message.textContent = `Pointer received: ${axis.name}, (${x.toFixed(5)}, ${y.toFixed(5)}). Updating; ${queue.length} queued.`;
        sendNext();
      });
      svg.append(rect);
    }
  };
  const acknowledge = content => {
    if (content.kind === 'ack') {
      pending = false; host.style.opacity = '1'; message.textContent = content.message;
      sendNext();
    }
  };
  model.on('change:svg', draw);
  model.on('msg:custom', acknowledge);
  draw();
  return () => {model.off('change:svg', draw); model.off('msg:custom', acknowledge);};
}};
