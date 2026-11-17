import io
import os
import subprocess
import tempfile
import zipfile
from pathlib import Path

from genbank_export import build_genbank_files


SNAPGENE_FORMATS = ["SnapGene", "SnapGene DNA"]


def _candidate_paths() -> list[Path]:
    paths: list[Path] = []
    env_path = (os.environ.get("SNAPGENE_EXE") or "").strip().strip('"')
    if env_path:
        paths.append(Path(env_path))
    paths.extend(
        [
            Path(r"C:\Program Files\SnapGene\SnapGene.exe"),
            Path(r"C:\Program Files\SnapGene.exe"),
            Path(r"C:\Program Files (x86)\SnapGene\SnapGene.exe"),
        ]
    )
    seen = set()
    unique = []
    for p in paths:
        key = str(p).lower()
        if key not in seen:
            unique.append(p)
            seen.add(key)
    return unique


def validate_snapgene_cli(path: Path) -> tuple[bool, str]:
    if not path.exists():
        return False, "文件不存在"
    try:
        proc = subprocess.run(
            [str(path), "--help"],
            capture_output=True,
            text=True,
            timeout=20,
            check=False,
        )
    except Exception as exc:
        return False, str(exc)
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode == 0 and (not output.strip() or ("convert" in output and "SnapGene" in output)):
        return True, "OK"
    return False, output.strip() or f"exit code {proc.returncode}"


def find_snapgene_cli() -> tuple[Path | None, str]:
    errors = []
    for path in _candidate_paths():
        ok, detail = validate_snapgene_cli(path)
        if ok:
            return path, detail
        errors.append(f"{path}: {detail}")
    return None, "；".join(errors)


def _convert_one(snapgene_exe: Path, gbk_name: str, gbk_data: bytes, out_dir: Path) -> tuple[str, bytes]:
    stem = Path(gbk_name).stem
    in_path = out_dir / gbk_name
    out_path = out_dir / f"{stem}.dna"
    in_path.write_bytes(gbk_data)

    last_error = ""
    for fmt in SNAPGENE_FORMATS:
        try:
            proc = subprocess.run(
                [str(snapgene_exe), "--accept-eula", "-q", "-c", fmt, "-i", str(in_path), "-o", str(out_path)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                timeout=120,
                check=False,
            )
        except subprocess.TimeoutExpired:
            last_error = "SnapGene CLI 转换超时，请确认 SnapGene 已激活且命令行转换不会弹出交互窗口。"
            continue
        if proc.returncode == 0 and out_path.exists() and out_path.stat().st_size > 0:
            return out_path.name, out_path.read_bytes()
        if proc.returncode == 0:
            last_error = "SnapGene CLI 正常退出但没有生成输出文件，请确认 SnapGene 已激活且当前版本支持命令行转换。"
        else:
            last_error = f"exit code {proc.returncode}"
        if out_path.exists():
            out_path.unlink()

    raise RuntimeError(f"SnapGene 转换失败：{last_error}")


def build_snapgene_dna_bundle(rows: list[dict]) -> tuple[str, bytes, str, Path]:
    snapgene_exe, detail = find_snapgene_cli()
    if snapgene_exe is None:
        raise RuntimeError(
            "当前环境未检测到 SnapGene Desktop CLI，已提供 GenBank `.gbk`；"
            "本地安装 SnapGene 后可导出 `.dna`。"
            f" 检测详情：{detail}"
        )

    gbk_files = build_genbank_files(rows)
    with tempfile.TemporaryDirectory(prefix="primer_snapgene_") as td:
        out_dir = Path(td)
        dna_files = [_convert_one(snapgene_exe, name, data, out_dir) for name, data in gbk_files]

    if len(dna_files) == 1:
        return dna_files[0][0], dna_files[0][1], "application/octet-stream", snapgene_exe

    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name, data in dna_files:
            zf.writestr(name, data)
    return "SnapGene_DNA_Maps.zip", bio.getvalue(), "application/zip", snapgene_exe
