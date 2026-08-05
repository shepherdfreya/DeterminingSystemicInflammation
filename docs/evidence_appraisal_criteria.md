## Evidence Appraisal Criteria

MR evidence was appraised using the following criteria.

### Instrument Selection

Exposure instruments were selected using:

- Genome-wide significance: `P < 5 x 10^-8`
- LD clumping: `R2 < 0.001` within a `10,000 kb` window
- Instrument strength: `F statistic > 10`

The F statistic was calculated as:

```text
F = (beta.exposure / se.exposure)^2
```

Proxy SNPs were searched for exposure instruments missing from the outcome GWAS using:

- Proxy LD threshold: `R2 > 0.8`
- Search window: `5,000 kb`
- Maximum proxy candidates: `5,000`

Where multiple proxies were available, the SNP with the highest `R2` was retained. Proxy SNPs were also required to pass `F statistic > 10`.

### Sensitivity Analyses

The following sensitivity analyses were included:

- Steiger filtering: SNPs were retained only where the inferred causal direction was exposure to outcome.
- Colocalization: evidence of a shared causal variant was assessed using `PP.H4.abf`.
- Heterogeneity: assessed using Cochran's Q.
- Directional pleiotropy: assessed using the MR-Egger intercept.

Colocalization was considered to pass where:

```text
PP.H4.abf > 0.8
```

Heterogeneity and pleiotropy tests were considered to pass where there was no strong evidence of violation:

```text
Cochran's Q p-value > 0.05
MR-Egger intercept p-value > 0.05
```

Where heterogeneity or pleiotropy tests were not applicable because too few SNPs were available, these were labelled as `Not applicable` and were not treated as failed criteria.

### Evidence Tiers

| Tier | pQTL Class | Bonferroni Significant | FDR Significant | Colocalization `H4 > 0.8` | Heterogeneity | Pleiotropy | Steiger Direction |
|---|---|---:|---:|---:|---:|---:|---:|
| Tier A | Cis | Yes | Yes | Yes | Pass or Not applicable | Pass or Not applicable | Yes |
| Tier B | Cis | No | Yes | Yes | Pass or Not applicable | Pass or Not applicable | Yes |
| Tier C | Cis/Trans | Not required | Yes | Yes | Pass or Not applicable | Pass or Not applicable | Yes |
| No Tier | Any | Criteria not met | Criteria not met | Criteria not met | Criteria not met | Criteria not met | Criteria not met |

Rows that did not meet the criteria for Tier A, Tier B, or Tier C were assigned `No Tier`.
