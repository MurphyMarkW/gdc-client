package main

import "os"

import cli "github.com/codegangsta/cli"
import log "github.com/Sirupsen/logrus"

import "github.com/gudCodes/gdc-client/client"
import "github.com/gudCodes/gdc-client/client/download"

// Deconstruct a cli context and call upload
// using context arguments and flags.
func upload(c *cli.Context) {
}

// Deconstruct a cli context and call query
// using context arguments and flags.
func query(c *cli.Context) {
}

// Logging flags shared between all commands.
var log_flags = []cli.Flag{
	cli.BoolFlag{
		Name:  "verbose",
		Usage: "verbose logging",
	},
	cli.BoolFlag{
		Name:  "debug",
		Usage: "debug logging",
	},
}

// Auth flags shared between all commands.
var auth_flags = []cli.Flag{
	cli.StringFlag{
		Name:  "T, token",
		Usage: "token string",
	},
	cli.StringFlag{
		Name:  "t, token-file",
		Usage: "token file",
	},
}

// Host flags shared between all commands.
var host_flags = []cli.Flag{
	cli.StringFlag{
		Name:   "H, host",
		Usage:  "GDC API host",
		Value:  "gdc-api.nci.nih.gov",
		EnvVar: "GDC_API_HOST",
	},
	cli.IntFlag{
		Name:   "P, port",
		Usage:  "GDC API port",
		Value:  443,
		EnvVar: "GDC_API_PORT",
	},
}

func setup(c *cli.Context) error {

	log.SetOutput(os.Stderr)

	log.SetLevel(log.WarnLevel)

	if c.Bool("verbose") {
		log.SetLevel(log.InfoLevel)
	}

	if c.Bool("debug") {
		log.SetLevel(log.DebugLevel)
	}

	return nil
}

func main() {
	app := cli.NewApp()

	app.Name = "gdc-client"
	app.Usage = "The Genomics Data Commons Command Line Client"
	app.Version = client.Version

	app.UsageText = "gdc-client [command] ..."

	var globals = []cli.Flag{}
	globals = append(globals, log_flags...)
	globals = append(globals, auth_flags...)
	globals = append(globals, host_flags...)

	app.Commands = []cli.Command{
		{
			Name:      "query",
			Flags:     append(globals, []cli.Flag{}...),
			Usage:     "query metadata using GDC DSL",
			Before:    setup,
			Action:    query,
			Category:  "metadata query",
			ArgsUsage: "QUERY",
		},
		{
			Name:      "upload",
			Flags:     append(globals, []cli.Flag{}...),
			Usage:     "upload data by GDC UUID",
			Before:    setup,
			Action:    upload,
			Category:  "data transfer",
			ArgsUsage: "ID [PATH]",
		},
		{
			Name:      "download",
			Flags:     append(globals, []cli.Flag{}...),
			Usage:     "download data to stdout",
			Before:    setup,
			Action:    download.DownloadAction,
			Category:  "data transfer",
			ArgsUsage: "ID",
		},
		{
			Name:      "download-bulk",
			Flags:     append(globals, download.Flags...),
			Usage:     "download data in bulk to filesystem",
			Before:    setup,
			Action:    download.DownloadBulkAction,
			Category:  "data transfer",
			ArgsUsage: "[ID]+",
		},
	}

	app.Run(os.Args)
}
