"""Real JupyterLab/Chromium pointer and widget tests, using a notebook copy.

Requires Playwright and a running localhost JupyterLab. See README_xsec_bin_helper.md.
Browser actions perform every edit. A second kernel channel only reads state and
checks counts; it does not invoke edit callbacks.
"""

from pathlib import Path
import argparse
import datetime
import json
import re
import time
import uuid
import requests
import websocket
from playwright.sync_api import sync_playwright, expect

ROOT = Path(__file__).resolve().parents[2]


class KernelReader:
    def __init__(self, base, token, notebook):
        sessions = requests.get(
            base + "/api/sessions",
            headers={"Authorization": "token " + token},
            timeout=15,
        ).json()
        kid = next(s["kernel"]["id"] for s in sessions if s["path"] == notebook)
        self.ws = websocket.create_connection(
            base.replace("http", "ws", 1)
            + f"/api/kernels/{kid}/channels?token={token}",
            origin=base,
            timeout=60,
        )

    def execute(self, code):
        ident = uuid.uuid4().hex
        self.ws.send(
            json.dumps(
                {
                    "header": {
                        "msg_id": ident,
                        "username": "validation",
                        "session": uuid.uuid4().hex,
                        "msg_type": "execute_request",
                        "version": "5.3",
                        "date": datetime.datetime.now(
                            datetime.timezone.utc
                        ).isoformat(),
                    },
                    "parent_header": {},
                    "metadata": {},
                    "channel": "shell",
                    "content": {
                        "code": code,
                        "silent": False,
                        "store_history": False,
                        "user_expressions": {},
                        "allow_stdin": False,
                        "stop_on_error": True,
                    },
                }
            )
        )
        output = ""
        while True:
            msg = json.loads(self.ws.recv())
            if msg.get("parent_header", {}).get("msg_id") != ident:
                continue
            if msg["msg_type"] == "error":
                raise AssertionError(msg["content"]["evalue"])
            if msg["msg_type"] == "stream":
                output += msg["content"]["text"]
            if (
                msg["msg_type"] == "status"
                and msg["content"]["execution_state"] == "idle"
            ):
                return output

    def state(self):
        return json.loads(self.execute("""print(json.dumps({
"vertices":dashboard.vertices,"counts":[len(f) for f in dashboard.selected_frames],
"revision":dashboard.revision,"valid":dashboard.valid,"geometry":dashboard.canvas.geometry,
"edges":{k:([e.tolist() for e in v] if k=="xB_by_Q2" else v.tolist()) for k,v in dashboard.edges.items()},
"limits":dashboard.limits,"settings":dashboard.settings}))"""))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default="http://127.0.0.1:8893")
    parser.add_argument("--token", required=True)
    parser.add_argument("--chromium", required=True)
    parser.add_argument("--case", choices=["both", "single"], default="both")
    args = parser.parse_args()
    import nbformat

    out = ROOT / "output/lt_bin_validation" / ("frontend_" + args.case)
    out.mkdir(parents=True, exist_ok=True)
    n = nbformat.read(ROOT / "src/xsec_extract/xsec_bin_helper.ipynb", as_version=4)
    if args.case == "single":
        n.cells[4].source = n.cells[4].source.replace(
            'KINEMATICS = ["KinC_x36_4", "KinC_x36_5_407"]',
            'KINEMATICS = ["KinC_x36_4"]',
        )
        n.cells[4].source = n.cells[4].source.replace(
            'REFERENCE_KINEMATIC = "KinC_x36_5_407"',
            'REFERENCE_KINEMATIC = "KinC_x36_4"',
        )
    # Collapsing sources improves navigation; all real notebook cells still run.
    for cell in n.cells:
        if cell.cell_type == "code":
            cell.metadata["jupyter"] = {"source_hidden": True}
    path = out / "interaction.ipynb"
    relative = str(path.relative_to(ROOT))
    # Only shut down this generated fixture's prior kernel, never a source notebook.
    auth = {"Authorization": "token " + args.token}
    for session in requests.get(
        args.url + "/api/sessions", headers=auth, timeout=15
    ).json():
        if session["path"] == relative:
            requests.delete(
                args.url + "/api/sessions/" + session["id"], headers=auth, timeout=15
            ).raise_for_status()
    nbformat.write(n, path)
    with sync_playwright() as p:
        browser = p.chromium.launch(
            executable_path=args.chromium,
            headless=True,
            args=["--no-sandbox", "--disable-dev-shm-usage"],
        )
        page = browser.new_page(
            viewport={"width": 1500, "height": 1400}, device_scale_factor=1
        )
        errors = []
        page.on("pageerror", lambda e: errors.append(str(e)))
        page.goto(
            args.url
            + "/lab/workspaces/lt-test-"
            + args.case
            + "/tree/"
            + relative
            + "?token="
            + args.token
        )
        page.get_by_text("Python 3 (ipykernel) | Idle", exact=False).wait_for(
            timeout=90000
        )
        page.get_by_text("Run", exact=True).first.click()
        run = page.locator(".lm-Menu-itemLabel").filter(
            has_text=re.compile("^Run All Cells$")
        )
        expect(run.locator("..")).not_to_have_class(
            re.compile("lm-mod-disabled"), timeout=60000
        )
        run.click()
        d = page.locator(".lt-dashboard")
        d.wait_for(timeout=240000)
        d.locator(".lt-canvas svg").wait_for(timeout=60000)
        reader = KernelReader(args.url, args.token, relative)
        reader.execute('print("kernel ready")')
        baseline = reader.state()
        results = []

        def snapshot(label):
            d.evaluate("el=>el.scrollIntoView({block:'center'})")
            page.wait_for_timeout(500)  # allow Jupyter's scroll/paint cycle to settle
            d.screenshot(path=str(out / (label + ".png")))
            page.screenshot(path=str(out / (label + "_jupyter.png")))

        def wait_revision(previous):
            for _ in range(80):
                state = reader.state()
                if state["revision"] > previous:
                    return state
                page.wait_for_timeout(100)
            raise AssertionError("Widget edit was not acknowledged by kernel")

        def click_corner(view, xb, q2, number):
            state = reader.state()
            axis = next(a for a in state["geometry"] if a["name"] == view)
            x, y = (
                (xb, q2)
                if view == "xB"
                else (q2, (0.9382720813**2 + q2 * (1 / xb - 1)) ** 0.5)
            )
            rx = (x - axis["xlim"][0]) / (axis["xlim"][1] - axis["xlim"][0])
            ry = (axis["ylim"][1] - y) / (axis["ylim"][1] - axis["ylim"][0])
            target = d.locator(f'[data-view="{view}"]')
            target.scroll_into_view_if_needed()
            box = target.bounding_box()
            target.click(position={"x": rx * box["width"], "y": ry * box["height"]})
            expect(d.locator(".lt-canvas [role=status]")).to_contain_text(
                f"Corner {number} accepted", timeout=30000
            )
            after = wait_revision(state["revision"])
            assert abs(after["vertices"][number - 1][0] - xb) < 0.002
            assert abs(after["vertices"][number - 1][1] - q2) < 0.02
            for imap in range(2):
                for i, v in enumerate(after["vertices"]):
                    if v is not None:
                        expect(d.locator(f'[id="corner-{imap}-{i+1}"]')).to_be_visible()
            results.append({"action": f"click {number} in {view}", "state": after})
            return after

        vertices = [(0.285, 3.25), (0.34, 3.3), (0.435, 4.7), (0.375, 4.6)]
        # Create a whole diamond from each 2D view, retaining all four labels.
        for view in ("xB", "W"):
            for i, v in enumerate(vertices):
                click_corner(view, *v, i + 1)
            snapshot(view + "_four_corners")
            original = reader.state()["vertices"]
            click_corner(view, 0.37, 4.5, 4)
            before = reader.state()["revision"]
            d.get_by_role("button", name="Undo", exact=True).click()
            restored = wait_revision(before)
            assert restored["vertices"] == original
            snapshot(view + "_undo_move4")
            before = restored["revision"]
            d.get_by_role("button", name="Clear", exact=True).click()
            cleared = wait_revision(before)
            assert (
                cleared["vertices"] == [None] * 4
                and cleared["counts"] == baseline["counts"]
            )
            expect(d.locator('[id^="corner-"]')).to_have_count(0)

        # Numeric fallback, arbitrary target corner, and controls through the DOM.
        state = reader.state()
        axis = next(a for a in state["geometry"] if a["name"] == "xB")
        target = d.locator('[data-view="xB"]')
        target.scroll_into_view_if_needed()
        box = target.bounding_box()
        for xb, q2 in vertices:
            x = box["x"] + box["width"] * (xb - axis["xlim"][0]) / (
                axis["xlim"][1] - axis["xlim"][0]
            )
            y = box["y"] + box["height"] * (axis["ylim"][1] - q2) / (
                axis["ylim"][1] - axis["ylim"][0]
            )
            page.mouse.click(x, y)
        expect(d.locator('[id="corner-1-4"]')).to_be_visible(timeout=30000)
        expect(d.locator(".lt-canvas [role=status]")).to_contain_text(
            "Corner 4 accepted", timeout=30000
        )
        state = reader.state()
        assert all(v is not None for v in state["vertices"])
        results.append({"action": "four rapid clicks queued", "state": state})
        snapshot("rapid_four_corners")
        before = state["revision"]
        d.get_by_role("button", name="Clear", exact=True).click()
        wait_revision(before)
        for i, v in enumerate(vertices):
            d.locator("select").first.select_option(label=str(i + 1))
            inputs = d.locator("input[type=number]")
            inputs.nth(0).fill(str(v[0]))
            inputs.nth(0).press("Tab")
            inputs.nth(1).fill(str(v[1]))
            inputs.nth(1).press("Tab")
            before = reader.state()["revision"]
            d.get_by_role("button", name="Set corner", exact=True).click()
            wait_revision(before)
        assert reader.state()["vertices"] == [list(v) for v in vertices]
        before = reader.state()["revision"]
        # Indices: corner x/Q2, three rows of min/max/bin count, phi/map/sparse.
        inputs = d.locator("input[type=number]")
        inputs.nth(11).fill("12")
        inputs.nth(11).press("Tab")
        state = wait_revision(before)
        assert len(state["edges"]["phi"]) == 13
        before = state["revision"]
        d.locator("select").nth(1).select_option("uniform")
        state = wait_revision(before)
        before = state["revision"]
        inputs.nth(4).fill("3")
        inputs.nth(4).press("Tab")
        state = wait_revision(before)
        assert len(state["edges"]["tprime"]) == 4
        # Click directly from the focused field: the button must capture its
        # visible text even if FloatText has not synced to the kernel yet.
        inputs.nth(2).fill("-.5")
        before = state["revision"]
        d.get_by_role("button", name="Apply outer limits", exact=True).click()
        state = wait_revision(before)
        assert state["limits"]["tprime"] == [-0.5, 0]
        expect(d.locator("[role=status]").last).to_contain_text("Outer limits applied")
        snapshot("changed_bins_limits")
        # Read-only independent physics audit after actual browser edits.
        reader.execute("""
for setting, f in zip(dashboard.settings, dashboard.selected_frames):
    raw=all_data[setting]
    mask=raw.mmiss_all.between(*MMISS_LIMITS)&raw.tprime.between(*TPRIME_LIMITS)&raw.Q2.between(*Q2_LIMITS)&raw.xB.between(*XB_LIMITS)
    polygon=dashboard._polygon()
    for a,b in zip(polygon,np.roll(polygon,-1,axis=0)):
        mask &= (b[0]-a[0])*(raw.Q2-a[1])-(b[1]-a[1])*(raw.xB-a[0])>=-1e-12
    assert set(raw.index[mask])==set(f.index)
    assert dashboard.occupancy_frame.query("setting == @setting").events.sum()==len(f)
assert len(data_selected)==sum(len(f) for f in dashboard.selected_frames)
""")
        reader.execute(f'dashboard.export_figure({str(out/"browser_selection")!r})')
        assert not errors, errors
        (out / "results.json").write_text(
            json.dumps(
                {
                    "case": args.case,
                    "browser": "Chromium",
                    "frontend": "JupyterLab",
                    "baseline": baseline["counts"],
                    "actions": results,
                    "final": state,
                    "javascript_errors": errors,
                    "checks": [
                        "four corners from each view",
                        "all labels on both maps",
                        "move corner 4",
                        "undo move",
                        "clear",
                        "rapid click queue",
                        "numeric corners",
                        "phi and tprime bin counts",
                        "edge mode",
                        "outer limits",
                        "independent mask/occupancy",
                    ],
                },
                indent=2,
            )
            + "\n"
        )
        reader.ws.close()
        browser.close()
        print("FRONTEND VALIDATION PASSED", out)


if __name__ == "__main__":
    main()
