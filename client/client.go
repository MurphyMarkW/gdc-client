package client

import "io"
import "fmt"
import "runtime"
import "net/http"
import "net/http/httputil"

import log "github.com/Sirupsen/logrus"
import errors "github.com/pkg/errors"

const AGENT = "gdc-client/%s (%s %s) %s (%s)"
const RESOURCE = "https://%s:%d%s"

type Client struct {
	http.Client
	host string
	port uint16
}

func NewClient(host string, port uint16) *Client {
	return &Client{
		host: host,
		port: port,
	}
}

func (c *Client) Request(verb string, path string, token string, body io.Reader) (*http.Request, error) {
	var url = fmt.Sprintf(RESOURCE, c.host, c.port, path)

	req, err := http.NewRequest(verb, url, body)
	if err != nil {
		return nil, errors.Wrap(err, "error generating new request")
	}

	agent := fmt.Sprintf(AGENT,
		Version,
		runtime.GOOS,
		runtime.GOARCH,
		runtime.Version(),
		runtime.Compiler,
	)

	req.Header.Set("User-Agent", agent)

	if len(token) > 0 {
		req.Header.Set("X-Auth-Token", token)
	}

	dump, err := httputil.DumpRequest(req, false)
	if err != nil {
		return nil, errors.Wrap(err, "error dumping HTTP request")
	}

	log.Debug(fmt.Sprintf("HTTP request:\n%s", dump))

	return req, err
}

func (c *Client) Do(req *http.Request) (*http.Response, error) {
	res, err := c.Client.Do(req)
	if err != nil {
		return nil, errors.Wrap(err, "error making HTTP request")
	}

	dump, err := httputil.DumpResponse(res, false)
	if err != nil {
		return nil, errors.Wrap(err, "error dumping HTTP response")
	}

	log.Debug(fmt.Sprintf("HTTP response:\n%s", dump))

	return res, err
}

func (c *Client) Head(path string, token string) (*http.Response, error) {
	req, err := c.Request(http.MethodHead, path, token, nil)
	if err != nil {
		return nil, errors.Wrap(err, "error generating HEAD request")
	}

	return c.Do(req)
}

func (c *Client) Get(path string, token string) (*http.Response, error) {
	req, err := c.Request(http.MethodGet, path, token, nil)
	if err != nil {
		return nil, errors.Wrap(err, "error generating GET request")
	}

	return c.Do(req)
}

func (c *Client) Put(path string, token string, body io.Reader) (*http.Response, error) {
	req, err := c.Request(http.MethodPut, path, token, body)
	if err != nil {
		return nil, errors.Wrap(err, "error generating PUT request")
	}

	return c.Do(req)
}

func (c *Client) Post(path string, token string, body io.Reader) (*http.Response, error) {
	req, err := c.Request(http.MethodPost, path, token, body)
	if err != nil {
		return nil, errors.Wrap(err, "error generating POST request")
	}

	return c.Do(req)
}

func (c *Client) Delete(path string, token string) (*http.Response, error) {
	req, err := c.Request(http.MethodDelete, path, token, nil)
	if err != nil {
		return nil, errors.Wrap(err, "error generating DELETE request")
	}

	return c.Do(req)
}
