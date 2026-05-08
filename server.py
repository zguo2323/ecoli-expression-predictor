import importlib.util
import inspect
import json
import sys
from pathlib import Path

from fastmcp import FastMCP

sys.path.insert(0, str(Path(__file__).parent))

mcp = FastMCP("BioE234 Expression Predictor")

MODULES_DIR = Path(__file__).parent / "modules"

TYPE_MAP = {"string": str, "integer": int, "number": float, "boolean": bool}


def _load_module(py_path: Path):
    spec = importlib.util.spec_from_file_location(py_path.stem, str(py_path))
    mod = importlib.util.module_from_spec(spec)
    sys.modules[py_path.stem] = mod
    spec.loader.exec_module(mod)
    return mod


def _find_mcp_class(mod):
    for name, obj in inspect.getmembers(mod, inspect.isclass):
        if name.endswith("MCP") and hasattr(obj, "initiate") and hasattr(obj, "run"):
            return obj
    return None


def _make_tool_fn(instance, inputs, tool_name, description):
    """Build a dynamic function with the right parameter signature for FastMCP."""
    param_names = [inp["name"] for inp in inputs]
    param_types = [TYPE_MAP.get(inp.get("type", "string"), str) for inp in inputs]

    def tool_fn(**kwargs):
        result = instance.run(**kwargs)
        if isinstance(result, (dict, list)):
            return json.dumps(result, indent=2)
        return str(result)

    sig_params = [
        inspect.Parameter(pname, inspect.Parameter.POSITIONAL_OR_KEYWORD, annotation=ptype)
        for pname, ptype in zip(param_names, param_types)
    ]
    tool_fn.__signature__ = inspect.Signature(sig_params)
    tool_fn.__annotations__ = {pname: ptype for pname, ptype in zip(param_names, param_types)}
    tool_fn.__name__ = tool_name
    tool_fn.__doc__ = description
    return tool_fn


registered = []

for module_dir in sorted(MODULES_DIR.iterdir()):
    if not module_dir.is_dir() or module_dir.name.startswith("_"):
        continue

    py_file = module_dir / f"{module_dir.name}.py"
    if not py_file.exists():
        continue

    mod = _load_module(py_file)
    mcp_class = _find_mcp_class(mod)
    if mcp_class is None:
        continue

    for json_file in sorted(module_dir.glob("*.json")):
        try:
            config = json.loads(json_file.read_text())
        except json.JSONDecodeError:
            continue

        exec_details = config.get("execution_details", {})
        mcp_name = exec_details.get("mcp_name")
        if not mcp_name:
            continue

        instance = mcp_class(config)
        instance.initiate()

        inputs = config.get("inputs", [])
        description = config.get("description", f"MCP tool: {mcp_name}")

        tool_fn = _make_tool_fn(instance, inputs, mcp_name, description)
        mcp.tool(name=mcp_name, description=description)(tool_fn)
        registered.append(mcp_name)


if __name__ == "__main__":
    print(f"[server] Registered {len(registered)} tools: {registered}", file=sys.stderr)
    mcp.run()
