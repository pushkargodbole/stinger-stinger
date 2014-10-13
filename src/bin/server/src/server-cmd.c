#include <stddef.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/un.h>

#include "stinger_core/xmalloc.h"

static char *
make_socket_name (const char *gname)
{
  char *dname, *fname;
  size_t dname_len, gname_len;

  if (!(dname = getenv ("TMPDIR")))
    if (!(dname = getenv ("TMP")))
      dname = "/tmp";

  dname_len = strlen (dname);
  gname_len = strlen (gname);
  fname = xmalloc ((dname_len + gname_len + 5) * sizeof (*fname));
  strcpy (fname, dname);
  strcpy (&fname[dname_len], gname);
  strcpy (&fname[dname_len+gname_len], "-ctl");
  return fname;
}

void
unlink_control_socket (const char *gname, int sock)
{
  char *fname = make_socket_name (gname);
  close (sock);
  unlink (fname);
  free (fname);
}

int
make_control_socket (const char *gname)
{
  char *fname;
  struct sockaddr_un name;
  int sock;
  size_t size;

  fname = make_socket_name (gname);

  sock = socket (PF_LOCAL, SOCK_STREAM, 0);
  if (sock < 0) {
    perror ("socket");
    exit (EXIT_FAILURE);
  }

  name.sun_family = AF_LOCAL;
  strncpy (name.sun_path, fname, sizeof (name.sun_path));
  name.sun_path[sizeof (name.sun_path) - 1] = '\0';
  size = SUN_LEN (&name);

  if (bind (sock, (struct sockaddr *) &name, size) < 0) {
    perror ("bind");
    exit (EXIT_FAILURE);
  }

  free (fname);

  return sock;
}

static int
connect_control_socket (const char *gname)
{
  char *fname;
  struct sockaddr_un name;
  int sock;
  size_t size;

  fname = make_socket_name (gname);

  sock = socket (PF_LOCAL, SOCK_STREAM, 0);
  if (sock < 0) {
    perror ("socket");
    exit (EXIT_FAILURE);
  }

  name.sun_family = AF_LOCAL;
  strncpy (name.sun_path, fname, sizeof (name.sun_path));
  name.sun_path[sizeof (name.sun_path) - 1] = '\0';
  size = SUN_LEN (&name);

  if (connect (sock, (struct sockaddr *) &name, size) < 0) {
    perror ("connect");
    exit (EXIT_FAILURE);
  }

  free (fname);

  return sock;
}

int
send_quit (const char *gname)
{
  int sock = connect_control_socket (gname);
  const char *cmd = "quit\n";
  ssize_t sz = write (sock, cmd, strlen (cmd));
  if (sz < 0) {
    perror ("control write");
    return -1;
  }
  close (sock);
  return 0;
}

static int quit;

static int
read_from_client (int fd)
{
  static char buf[1024];
  ssize_t len;

  len = read (fd, buf, 1000);
  if (len < 0) {
    perror ("client read");
  } else if (len > 0) {
    /* fprintf (stderr, "CLIENT COMMAND \"%s\"\n", buf); */
    /* Handle the only command... */
    if (buf[0] == 'q')
      quit = 1;
  }
  return len;
}

int
server_loop (int sock)
{
  fd_set active_fd_set, read_fd_set;

  quit = 0;

  if (listen (sock, 5) < 0) {
    perror ("listen");
    return -1;
  }

  FD_ZERO (&active_fd_set);
  FD_SET (sock, &active_fd_set);

  while (!quit) {
    read_fd_set = active_fd_set;
    if (select (FD_SETSIZE, &read_fd_set, NULL, NULL, NULL) < 0) {
      perror ("select");
      return -1;
    }

    for (int i = 0; i < FD_SETSIZE; ++i)
      if (FD_ISSET (i, &read_fd_set)) {
        if (i == sock) {
          struct sockaddr_un client;
          socklen_t size;
          int new;
          size = sizeof (client);
          new = accept (sock, (struct sockaddr *) &client, &size);
          if (new < 0) {
            perror ("accept");
            return -1;
          }
          FD_SET (new, &active_fd_set);
        }
        else {
          if (read_from_client (i) <= 0) {
            close (i);
            FD_CLR (i, &active_fd_set);
          }
        }
      }
  }
}
