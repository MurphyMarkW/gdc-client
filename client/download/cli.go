package download

import "io"
import "os"
import "fmt"
import "mime"
import "path"
import "crypto"
import "crypto/md5"

import bar "gopkg.in/cheggaaa/pb.v1"
import cli "github.com/urfave/cli"
import log "github.com/Sirupsen/logrus"

var Flags = []cli.Flag{
	cli.StringFlag{
		Name:  "p, path",
		Usage: "path to target directory",
		Value: ".",
	},
	cli.StringFlag{
		Name:  "m, manifest",
		Usage: "GDC manifest file",
	},
}

// Verifies the contents of an io.Reader in a streaming fashion.
// If a hashing error is encountered when the proxied reader is
// closed, then an error is raised instead.
func verify(reader io.Reader, hash crypto.Hash, digest []byte) io.Reader {
	// TODO IMPLEMENT ME
	return reader
}

// Prints a progress meter for an io.Reader as it is consumed.
func progress(reader io.Reader, length int64) io.Reader {
	// TODO IMPLEMENT ME
	return reader
}

func DownloadAction(c *cli.Context) {
	if len(c.Args()) > 1 {
		log.Fatalf("more than one ID provided")
	}

	host := c.String("host")
	port := uint16(c.Int("port"))
	uuid := c.Args()[0]
	token := c.String("token")

	log.Info(fmt.Sprintf("Downloading %s...", uuid))

	client := NewDownloadClient(host, port)
	res, err := client.Download(uuid, token)
	if err != nil {
		log.Errorf("%s", err)
		return
	}

	io.Copy(os.Stdout, res.Body)
}

func DownloadBulkAction(c *cli.Context) {
	host := c.String("host")
	port := uint16(c.Int("port"))
	uuids := c.Args()
	token := c.String("token")

	if len(uuids) < 1 {
		log.Fatalf("must provide at least one ID")
	}

	log.Info("Downloading bulk...")
	for _, uuid := range uuids {
		log.Info(uuid)

		var result []byte

		client := NewDownloadClient(host, port)
		res, err := client.Download(uuid, token)
		if err != nil {
			log.Errorf("%s", err)
			return
		}

		content_disposition := res.Header.Get("Content-Disposition")
		media, params, err := mime.ParseMediaType(content_disposition)
		if err != nil {
			log.Errorf("%s", err)
			return
		}

		if media != "attachment" {
			log.Fatalf("non-attachment response from api")
		}

		if err := os.MkdirAll(uuid, 0700); err != nil {
			log.Fatalf(fmt.Sprintf("could not create directory %s", uuid))
		}

		ofs, err := os.Create(path.Join(uuid, params["filename"]))
		if err != nil {
			log.Errorf("%s", err)
			return
		}

		defer ofs.Close()

		// TODO what follows is a mess of progress metering and
		// hash generation / verification and needs cleaning up

		pbar := bar.New64(res.ContentLength).SetUnits(bar.U_BYTES)

		pbar.Prefix(uuid)
		pbar.Output = os.Stderr

		pbar.Start()

		proxy := pbar.NewProxyReader(res.Body)

		hash := md5.New()

		verified := io.TeeReader(proxy, hash)

		if _, err := io.Copy(ofs, verified); err != nil {
			log.Fatalf("%v", err)
		}

		pbar.Finish()

		if res.Header.Get("Content-MD5") != fmt.Sprintf("%x", hash.Sum(result)) {
			log.Errorf(fmt.Sprintf("received %x", hash.Sum(result)))
			log.Errorf(fmt.Sprintf("expected %s", res.Header.Get("Content-MD5")))
			log.Fatalf(fmt.Sprintf("error while hashing %s", uuid))
		}
	}
}
