# The GDC Client - Community Go Edition!
[![GitHub release](https://img.shields.io/github/release/gudcodes/gdc-client.svg)](https://github.com/MurphyMarkW/gdc-client/releases) [![MIT license](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE) [![Build Status](https://travis-ci.org/MurphyMarkW/gdc-client.svg?branch=master)](https://travis-ci.org/MurphyMarkW/gdc-client)

A community go implementation of the [Genomic Data Commons](https://gdc.nci.nih.gov/) Command Line Client and related libraries.

## Installation

Pre-built versions of the GDC Client [are available](https://github.com/MurphyMarkW/gdc-client/releases)
for most major operating systems and architectures. Simply download and
extract the binary in order to use the client.

### Building From Source

The GDC Client can be built from source using the Golang [Go command](https://golang.org/doc/).

```
go get github.com/MurphyMarkW/gdc-client
```

This will download and build the `gdc-client` package along with all external
dependencies.

## Usage

The Community GDC Client provides a collection of commands grouped by function.
A list of commands available and how to use them is available by running the
client without any arguments, or by specifying the `-h` or `--help` flags.

```
$ gdc-client
NAME:
   gdc-client - The Genomics Data Commons Command Line Client

USAGE:
   gdc-client [command] ...

VERSION:
   0.0.3

COMMANDS:
     help, h  Shows a list of commands or help for one command
   data transfer:
     upload         upload data by GDC UUID
     download       download data to stdout
     download-bulk  download data in bulk to filesystem
   metadata query:
     query  query metadata using GDC DSL

GLOBAL OPTIONS:
   --help, -h     show help
   --version, -v  print the version
```

The version of the client can be retrieved by supplying the `-v` or `--version`
flags without any command specified.

```
$ gdc-client --version
gdc-client version 0.0.3
```

Help for specific commands can be found by specifying the same `-h` or `--help`
flags after the command's name.

```
$ gdc-client download -h
NAME:
   gdc-client download - download data to stdout

USAGE:
   gdc-client download [command options] ID

CATEGORY:
   data transfer

OPTIONS:
   --verbose                     verbose logging
   --debug                       debug logging
   -T value, --token value       token string
   -t value, --token-file value  token file
   -H value, --host value        GDC API host (default: "gdc-api.nci.nih.gov") [$GDC_API_HOST]
   -P value, --port value        GDC API port (default: 443) [$GDC_API_PORT]
```

### Encryption

The GDC Client is shipped with SSL / TLS encryption built-in and enabled.
All pre-built versions of the GDC Client enforce verification of GDC API
certificates and use of SSL / TLS encrypted connections.

**Note** that proper SSL / TLS encryption and certificate verification can
be disabled when modifying and building the client from source. Doing so may
violate government regulations regarding the proper handling of protected data.
Please contact GDC support if you are unsure whether a modification to the
client is compliant with government regulations applicable to use of the GDC.

### Authentication & Authorization

The GDC API uses a token-based authentication mechanism for identifying users.
Access to protected data via the GDC API will require a GDC auth token. This
token identifies you as an individual and allows the GDC to track all actions
taken by the holder of this token. A GDC auth token can be acquired by visiting
the [GDC Portal](https://gdc-portal.nci.nih.gov/), logging in, and downloading
a token file from your account drop-down menu.

**Please secure your token and do not share with any other individual.**

A GDC token is meant to be tied to a single individual. Sharing your token
increases the likelihood of it being compromised and used by an unauthorized
party. In the event that a token is believed to have been compromised, please
contact [GDC Support](mailto:support@nci-gdc.datacommons.io) for assistance
in revoking the token and getting a report of data access for your review.

### Global Flags

The GDC Client operates as a thin HTTP(S) wrapper with some additional logic
layered on top. All commands interact with the GDC API in some fashion, so
certain flags are available for all commands. These include:

###### --verbose

Enables verbose logging. This can include useful information when used in
an interactive fashion.

###### --debug

Enables debug logging. This implies `--verbose` and includes information used
for identifying and diagnosing bugs in the client. Note that this can have a
negative impact on performance and may print sensitive information in plain
text. *Caution should be taken when used in production environments.*

###### -H / --host

In some cases it may be necessary to direct the client towards a different
host than the public GDC API. This may include certain types of HTTP(S)
proxies, such as UDT-enabled [Parcel](https://github.com/LabAdvComp/parcel)
proxies for long-distance, high latency networks.

###### -P / --port

As with the `-H` / `--host` flags, alternative ports other than the standard
`443` used for SSL-encrypted HTTP may be necessary.

###### -t / --token-file

A token file containing a GDC authentication token to be used with requests.
Protected data stored by the GDC will need an authentication token in order
to be accessed by clients. Authentication tokens can be downloaded from the
GDC portal after logging in.

### Environment Variables

The GDC Client is aware of certain environment variables.

###### GDC\_API\_HOST

The `GDC_API_HOST` environment variable performs the same action as
`-H` / `--host`. In the event that both `GDC_API_HOST` and `-H` / `--host`
are specified, the `-H` / `--host` will be used.

###### GDC\_API\_PORT

The `GDC_API_PORT` environment variable performs the same action as
`-P` / `--port`. In the event that both `GDC_API_PORT` and `-P` / `--port`
are specified, the `-P` / `--port` will be used.

###### HTTP{,S}\_PROXY / http{,s}\_proxy

The `HTTP_PROXY` and `HTTPS_PROXY` environment variables are both respected
by the GDC Client. Note that `HTTPS_PROXY` takes precedence over `HTTP_PROXY`
in all cases.

###### NO\_PROXY / no\_proxy

The `NO_PROXY` environment variable is respected by the GDC Client. This can
be used to avoid proxying connections to `gdc-api.nci.nih.gov` by adding the
GDC API hostname to the `NO_PROXY` environment variable:

```
$ export NO_PROXY=gdc-api.nci.nih.gov,${NO_PROXY}
```

This can be useful in environments that require low-volume traffic to be
proxied for security purposes, but allow certain high-volume traffic to known,
trusted hosts to bypass the proxy.

### Download

The `download` command is the most fundamental of the client functionality.
By providing the client with a GDC UUID, the client will make a secure request
to the GDC API and stream the data for the specified GDC UUID to standard out.

```
$ gdc-client download -h
NAME:
   gdc-client download - download data by GDC UUID

USAGE:
   gdc-client download [command options] id [path]

CATEGORY:
   data transfer

OPTIONS:
   --verbose                        verbose logging
   --debug                          debug logging
   -T, --token                      token string
   -t, --token-file                 token file
   -H, --host "gdc-api.nci.nih.gov" GDC API host [$GDC_API_HOST]
   -P, --port "443"                 GDC API port [$GDC_API_PORT]
```

**NOTE** that this command can produce significant output, especially when
working with large objects. This is an intended feature of the client that
facilitates the use of the client in downstream workflows. As an example,
the header from any BAM file in the GDC could be retrieved by piping the
output of `download` to `samtools`:

```
$ gdc-client download -t SUPER_SECRET_TOKEN_FILE e045df57-b8c4-4687-97bb-63b9f7a0357b | samtools view -H
@HD VN:1.5  SO:coordinate
@SQ SN:chr1 LN:248956422
@SQ SN:chr2 LN:242193529
@SQ SN:chr3 LN:198295559
@SQ SN:chr4 LN:190214555
@SQ SN:chr5 LN:181538259
@SQ SN:chr6 LN:170805979
@SQ SN:chr7 LN:159345973
@SQ SN:chr8 LN:145138636
@SQ SN:chr9 LN:138394717
...
```

### Download (Bulk)

The `download-bulk` command is a convenience function for the download of
multiple files simultaneously from the GDC. By providing one or more GDC
UUIDs, the client will make secure requests to the GDC API and stream the
data for the specified GDC UUIDs to the local filesystem. Directories and
files will be created based upon the GDC UUIDs and filenames reported by
the GDC. An optional target `path` may be supplied to specify the desired
root directory to be used for download.

```
$ gdc-client download-bulk -h
NAME:
   gdc-client download-bulk - download data in bulk to filesystem

USAGE:
   gdc-client download-bulk [command options] [ID]+

CATEGORY:
   data transfer

OPTIONS:
   --verbose                        verbose logging
   --debug                          debug logging
   -T, --token                      token string
   -t, --token-file                 token file
   -H, --host "gdc-api.nci.nih.gov" GDC API host [$GDC_API_HOST]
   -P, --port "443"                 GDC API port [$GDC_API_PORT]
   -p, --path "."                   path to target directory
   -m, --manifest                   GDC manifest file
```

**NOTE** this command assumes a POSIX-like file system on UNIX-like operating
systems, including OS X. On windows operating systems, the current or target
folder must support folder creation and file write operations for the current
user.

```
$ gdc-client download-bulk --verbose -t SUPER_SECRET_TOKEN_FILE e045df57-b8c4-4687-97bb-63b9f7a0357b
INFO[0000] Downloading bulk...
INFO[0000] e045df57-b8c4-4687-97bb-63b9f7a0357b
e045df57-b8c4-4687-97bb-63b9f7a0357b 22.69 MB / 22.69 MB [=============================] 100.00% 20s

$ ls -F
e045df57-b8c4-4687-97bb-63b9f7a0357b/

$ ls -F e045df57-b8c4-4687-97bb-63b9f7a0357b
113075.bam

$ samtools view -H e045df57-b8c4-4687-97bb-63b9f7a0357b/113075.bam
@HD VN:1.5  SO:coordinate
@SQ SN:chr1 LN:248956422
@SQ SN:chr2 LN:242193529
@SQ SN:chr3 LN:198295559
@SQ SN:chr4 LN:190214555
@SQ SN:chr5 LN:181538259
@SQ SN:chr6 LN:170805979
@SQ SN:chr7 LN:159345973
@SQ SN:chr8 LN:145138636
@SQ SN:chr9 LN:138394717
...

$ samtools index e045df57-b8c4-4687-97bb-63b9f7a0357b/113075.bam

$ samtools view e045df57-b8c4-4687-97bb-63b9f7a0357b/113075.bam chr1:412892-412892
SOLEXA7_0100:5:100:5593:4686    16  chr1    412888  0   20M *   0   0   ACTGACTGACTGACTGACTG    9.(9B7/*?B'@-=+03(>,
```
