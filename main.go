package main

import "os"

import "github.com/codegangsta/cli"
import "github.com/NCI-GDC/gdc-client/client/actions"

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

func main() {
	app := cli.NewApp()

	app.Name = "gdc-client"
	app.Usage = "The Genomics Data Commons Command Line Client"
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
			Action:    query,
			Category:  "metadata query",
			ArgsUsage: "query",
		},
		{
			Name:      "upload",
			Flags:     append(globals, []cli.Flag{}...),
			Usage:     "upload data by GDC UUID",
			Action:    upload,
			Category:  "data transfer",
			ArgsUsage: "id [path]",
		},
		{
			Name:      "download",
			Flags:     append(globals, []cli.Flag{}...),
			Usage:     "download data by GDC UUID",
			Action:    download,
			Category:  "data transfer",
			ArgsUsage: "id [path]",
		},
	}

	app.Run(os.Args)
}

// Deconstruct a cli context and call download
// using context arguments and flags.
func download(c *cli.Context) {
	actions.Download(
		c.String("host"),
		c.Int("port"),
		c.String("token"),
		c.Args()[0],
		c.Args()[1],
	)
}

// Deconstruct a cli context and call upload
// using context arguments and flags.
func upload(c *cli.Context) {
	actions.Upload(
		c.String("host"),
		c.Int("port"),
		c.String("token"),
		c.Args()[0],
		c.Args()[1],
	)
}

// Deconstruct a cli context and call query
// using context arguments and flags.
func query(c *cli.Context) {
	actions.Query(
		c.String("host"),
		c.Int("port"),
		c.String("token"),
		c.Args()[0],
	)
}
