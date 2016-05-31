The Genomic Data Commons Command Line Client Specifications
===

# DESCRIPTION

The GDC Client is a set of client side libraries and standalone scripts that, combined, allow for user interaction with the GDC API.

# REQUIREMENTS

The client must provide a light weight means of securely making arbitrary HTTPS calls to the GDC API.
All such calls must implement GDC security protocols.

The client must provide a library implementing common interactions with the GDC API including, but not limited to:
    - download
    - upload
    - search

The client must provide a command line client implementing common interactions with the GDC API including, but not limited to:
    - download
    - upload
    - search

# SECURITY

All calls made by the client must satisfy the following security requirements:
    - all traffic must be SSL / TLS encrypted with an approved cipher
    - all requests must verify the GDC certificate presented by the API
    - when provided by the user, all requests must provide the user token in form of an X-Auth-Header request header

## FEATURE - DOWNLOAD

Download should take as input one or more GDC UUIDs and produce one or more files.

Download must support resumption of partially completed downloads when downloading to disk.
Download must support MD5 checksum verification when downloading to disk.
Download should support common data-interchange formats as input, including but not limited to:
    - JSON
    - YAML

## FEATURE - QUERY

Queries should take as input one or more GDC-formatted queries and produce a JSON response.
