# The GDC Client - Community Go Edition!
[![Version](https://img.shields.io/badge/version-2016.04.0-red.svg?style=flat)](README.md) [![MIT license](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

A community go implementation of the Genomic Data Commons Command Line Client and related libraries.

References :

* https://gdc.nci.nih.gov/

## Installation

```
go get github.com/gudCodes/gdc-client
```

Dependencies :
* [github.com/codegangsta/cli](https://github.com/codegangsta/cli)

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
   0.0.1

COMMANDS:
  data transfer:
    upload      upload data by GDC UUID
    download        download data to stdout
    download-bulk   download data in bulk to filesystem

  metadata query:
    query   query metadata using GDC DSL

GLOBAL OPTIONS:
   --help, -h       show help
   --version, -v    print the version
```

Help for specific commands can be found by specifying the same `-h` or `--help`
flags after the command's name.

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

### Download

The `download` command is the most fundamental of the client functionality.
By providing the client with a GDC UUID, the client will make a secure request
to the GDC API and stream the data for the specified GDC UUID to standard out.

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

### BAM

**NOT IMPLEMENTED**

Per byte, the GDC stores more BAM formatted data than any other data type.
Because of this, certain functionality to ease the retrieval of regions of
interest from BAM files is included with the client. By paring download with
ranged requests, this can be used to slice small regions from otherwise large
BAM files:

```
$ gdc-client bam -t SUPER_SECRET_TOKEN_FILE e045df57-b8c4-4687-97bb-63b9f7a0357b chr1:412892-412892 | samtools view
SOLEXA7_0100:5:100:5593:4686    16  chr1    412888  0   20M *   0   0   ACTGACTGACTGACTGACTG    9.(9B7/*?B'@-=+03(>,
```
