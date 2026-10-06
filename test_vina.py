from docking.vina_runner import run_vina


box = {
    "center_x": 328.609,
    "center_y": 213.901,
    "center_z": 192.315,
    "size_x": 24.614,
    "size_y": 28.684,
    "size_z": 29.093
}


run_vina(
    ".\\vina.exe",
    "results/test_CDK1_receptor.pdbqt",
    "results/test_flavopiridol.pdbqt",
    "results/test_docking.pdbqt",
    box,
    exhaustiveness=16,
    num_modes=10
)