import os
import pandas as pd


def _strip_quotes(val):
    if isinstance(val, str):
        return val.strip('"').strip()
    return val


def load_promoter_table(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    df["Promoter"] = df["Promoter"].apply(_strip_quotes)
    df["Sequence"] = df["Sequence"].apply(_strip_quotes)
    df = df.rename(columns={
        "Promoter": "promoter_id",
        "Sequence": "sequence",
        "mean.RNA": "mean_RNA",
        "mean.prot": "mean_prot",
        "TSS.best": "TSS_best",
    })
    df = df.dropna(subset=["sequence"])
    return df[["promoter_id", "sequence", "mean_RNA", "mean_prot", "TSS_best"]]


def load_rbs_table(path: str) -> pd.DataFrame:
    df = pd.read_excel(path)
    df["RBS"] = df["RBS"].apply(_strip_quotes)
    df["Sequence"] = df["Sequence"].apply(_strip_quotes)
    df = df.rename(columns={
        "RBS": "rbs_id",
        "Sequence": "sequence",
        "mean.RNA": "mean_RNA",
        "mean.prot": "mean_prot",
        "mean.xlat": "mean_xlat",
    })
    df = df.dropna(subset=["sequence"])
    return df[["rbs_id", "sequence", "mean_RNA", "mean_prot", "mean_xlat"]]


def load_construct_table(path: str) -> pd.DataFrame:
    df = pd.read_excel(path, sheet_name="Constructs")
    df["Promoter"] = df["Promoter"].apply(_strip_quotes)
    df["RBS"] = df["RBS"].apply(_strip_quotes)
    df = df[(df["bad.prot"] == False) & (df["bad.RNA"] == False)]
    df = df.rename(columns={
        "Promoter": "promoter_id",
        "RBS": "rbs_id",
    })
    df = df.dropna(subset=["prot", "RNA"])
    return df[["promoter_id", "rbs_id", "RNA", "prot", "deltaG"]]


def build_dataset(
    sd01_path: str,
    sd02_path: str,
    sd03_path: str,
    output_path: str = "data/processed/constructs.parquet",
) -> pd.DataFrame:
    promoters_df = load_promoter_table(sd01_path)[["promoter_id", "sequence"]].rename(
        columns={"sequence": "promo_seq"}
    )
    rbs_df = load_rbs_table(sd02_path)[["rbs_id", "sequence"]].rename(
        columns={"sequence": "rbs_seq"}
    )
    constructs_df = load_construct_table(sd03_path)

    df = constructs_df.merge(promoters_df, on="promoter_id", how="left")
    df = df.merge(rbs_df, on="rbs_id", how="left")
    df = df.dropna(subset=["promo_seq", "rbs_seq"])

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_parquet(output_path, index=False)
    return df


class PipelineMCP:
    def __init__(self, config):
        self.config = config

    def initiate(self):
        pass

    def run(self, **kwargs):
        mcp_name = self.config.get("execution_details", {}).get("mcp_name")
        if mcp_name == "build_dataset":
            result = build_dataset(
                sd01_path=kwargs["sd01_path"],
                sd02_path=kwargs["sd02_path"],
                sd03_path=kwargs["sd03_path"],
            )
            return {"output_path": "data/processed/constructs.parquet"}
        raise ValueError(f"Unknown mcp_name: {mcp_name}")
