#if !defined(SERVER_CMD_HEADER_)
#define SERVER_CMD_HEADER_

int make_control_socket (const char *);
void unlink_control_socket (const char *, int);
int send_quit (const char *);
int server_loop (int);

#endif /* SERVER_CMD_HEADER_ */
