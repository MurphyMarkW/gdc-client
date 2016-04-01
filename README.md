# The GDC Client - Go Edition!
[![Version](https://img.shields.io/badge/version-2016.04.0-red.svg?style=flat)](README.md) [![MIT license](http://img.shields.io/badge/license-MIT-blue.svg?style=flat)](LICENSE)

Go implementation of the Genomic Data Commons Command Line Client and related libraries.

References :

* https://gdc.nci.nih.gov/

## Installation

```
go get github.com/gudCodes/gdc-client
```

Dependencies :
* [github.com/codegangsta/cli](https://github.com/codegangsta/cli)

## Usage

```
$ gdc-client
NAME:
   gdc-client - The Genomics Data Commons Command Line Client

USAGE:
   gdc-client [command] ...

VERSION:
   0.0.0

COMMANDS:
  data transfer:
    upload      upload data by GDC UUID
    download    download data by GDC UUID

  metadata query:
    query       query metadata using GDC DSL

GLOBAL OPTIONS:
   --help, -h       show help
   --version, -v    print the version
```

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
