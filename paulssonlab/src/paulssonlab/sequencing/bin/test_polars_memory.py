import polars as pl

dfs = []
for i in range(100):
    df = pl.DataFrame([{"string": "ATCG" * 2_000}] * 100_000).lazy()
    # df = pl.DataFrame([{"string": [12]*2_00}]*100_000).lazy()
    df = df.head(1000)
    df = df.collect()
    print("COLLECTED", i)
    dfs.append(df)
    # dfs.append(pl.from_arrow(df.to_arrow()))
# large_df = pl.concat(dfs)
# large_df = large_df.collect()
