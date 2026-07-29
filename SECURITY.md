# Security Policy

## Supported versions

The `0.1.x` line is experimental. Only the latest published `0.1.x` patch and
the default development branch receive security fixes; earlier `0.1.x` patches
may be superseded rather than patched in place. Versions older than `0.1.0` are
unsupported.

The API may change in later `0.x` minor releases, while `0.1.x` patch releases
may correct numerical behavior. A scientific-correctness issue is not
automatically a security vulnerability, but an incorrect certificate or an
exploitable resource-exhaustion path should be reported privately.

## Reporting a vulnerability

Send a private report to `jmartin@uliege.be` with:

- affected commit or version;
- reproducible steps or a minimal example;
- expected impact and threat model;
- any suggested mitigation;
- whether and when you intend to disclose publicly.

Do not open a public issue for an unpatched vulnerability or include secrets,
personal data, proprietary solver files, or malicious payloads beyond what is
needed to reproduce the problem. If encrypted reporting becomes available, this
file will list the key and channel.

The maintainer will acknowledge and assess reports as capacity permits, avoid
public disclosure before a reasonable remediation opportunity, and credit the
reporter if desired. No fixed response-time guarantee is made at this stage.

## Scientific correctness

Incorrect certification, silently ignored solver failure, unexpected
densification, denial of service from adversarial dimensions, unsafe
deserialization, and global RNG side effects can all be security or integrity
concerns in scientific workflows. Report them privately when exploitation or
confidentiality is involved; ordinary numerical bugs may use the public issue
tracker once one exists.
