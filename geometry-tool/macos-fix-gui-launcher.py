from pathlib import Path
import shutil
import sys

exe = shutil.which("hlb-gmy-gui")
with open(exe) as f:
    text = f.read()

shebang, body = text.split("\n", 1)
assert shebang[:2] == "#!"
script_interp = Path(shebang[2:])

current_interp = Path(sys.executable)
assert script_interp.resolve() == current_interp.resolve()

gui_interp = shutil.which("pythonw")
if gui_interp is None:
    # Try the conda python.app
    conda = os.environ["CONDA_PREFIX"]
    gui_interp = Path(conda) / "python.app/Contents/MacOS/python"
    assert gui_interp.exists()

new_text = f"""#!{gui_interp}
{body}"""

with open(exe, "w") as f:
    f.write(new_text)
