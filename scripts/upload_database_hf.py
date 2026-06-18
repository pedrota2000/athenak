from huggingface_hub import login, HfApi


login()


api = HfApi()
api.upload_large_folder(
    folder_path="/mnt/ceph/users/ptarancon/runs/turb_128_run_20260223_130640/vtk",
    repo_id="pedrota2000/NS_simulation",
    repo_type="dataset",
)