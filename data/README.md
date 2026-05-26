# Data

Place the whistle sub-clustering CSV at:

```
data/AllWhistlesSubClustering_final.csv
```

## Required columns

| Column | Description |
|--------|-------------|
| `recording` | Recording identifier |
| `whistle_type_chr` | Whistle sub-category label |
| `start_time` | Whistle onset (seconds) |
| `end_time` | Whistle offset (seconds) |
| `recording_duration` | Duration of the recording (seconds) |
| `whistle_type` | Numeric whistle type ID |
| `whistle_name` | Whistle name |

Rows with `Unknown` in `whistle_type_chr` are excluded by the analysis scripts.
