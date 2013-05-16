#include <stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/wait.h>
#include <sys/types.h>

static char **cmd2argv(const char *cmd)
{
	int i, beg, end, argc;
	char **argv, *str;
	end = strlen(cmd);
	for (i = end - 1; i >= 0; --i)
		if (!isspace(cmd[i])) break;
	end = i + 1;
	for (beg = 0; beg < end; ++beg)
		if (!isspace(cmd[beg])) break;
	if (beg == end) return 0;
	for (i = beg + 1, argc = 0; i < end; ++i)
		if (isspace(cmd[i]) && !isspace(cmd[i-1]))
			++argc;
	argv = (char**)calloc(argc + 2, sizeof(void*));
	argv[0] = str = (char*)calloc(end - beg + 1, 1);
	strncpy(argv[0], cmd + beg, end - beg);
	for (i = argc = 1; i < end - beg; ++i)
		if (isspace(str[i])) str[i] = 0;
		else if (str[i] && str[i-1] == 0) argv[argc++] = &str[i];
	return argv;
}

#define KO_STDIN    1
#define KO_FILE     2
#define KO_PIPE     3
#define KO_HTTP     4
#define KO_FTP      5

typedef struct {
	int type, fd;
	pid_t pid;
} koaux_t;

void *kopen(const char *fn, int *_fd)
{
	koaux_t *aux = 0;
	*_fd = -1;
	if (strcmp(fn, "-") == 0) {
		aux = calloc(1, sizeof(koaux_t));
		aux->type = KO_STDIN;
		aux->fd = STDIN_FILENO;
	} else {
		const char *p, *q;
		for (p = fn; *p; ++p)
			if (!isspace(*p)) break;
		if (*p == '<') { // pipe open
			int need_shell, pfd[2];
			pid_t pid;
			// a simple check to see if we need to invoke a shell; not always working
			for (q = p + 1; *q; ++q)
				if (ispunct(*q) && *q != '.' && *q != '_' && *q != '-' && *q != ':')
					break;
			need_shell = (*q != 0);
			if (pipe(pfd) != 0) return 0;
			pid = vfork();
			if (pid == -1) { /* vfork() error */
				close(pfd[0]); close(pfd[1]);
				return 0;
			}
			if (pid == 0) { /* the child process */
				char **argv; /* FIXME: I do not know if this will lead to a memory leak */
				close(pfd[0]);
				dup2(pfd[1], STDOUT_FILENO);
				close(pfd[1]);
				if (!need_shell) {
					argv = cmd2argv(p + 1);
					execvp(argv[0], argv);
					free(argv[0]); free(argv);
				} else execl("/bin/sh", "sh", "-c", p + 1, NULL);
				exit(1);
			} else { /* parent process */
				close(pfd[1]);
				aux = calloc(1, sizeof(koaux_t));
				aux->type = KO_PIPE;
				aux->fd = pfd[0];
				aux->pid = pid;
			}
		} else {
#ifdef _WIN32
			*_fd = open(fn, O_RDONLY | O_BINARY);
#else
			*_fd = open(fn, O_RDONLY);
#endif
			if (*_fd >= 0) {
				aux = calloc(1, sizeof(koaux_t));
				aux->type = KO_FILE;
				aux->fd = *_fd;
			}
		}
	}
	if (aux) *_fd = aux->fd;
	return aux;
}

int kclose(void *a)
{
	koaux_t *aux = (koaux_t*)a;
	if (aux->type == KO_PIPE) {
		int status;
		pid_t pid;
		pid = waitpid(aux->pid, &status, WNOHANG);
		if (pid != aux->pid) kill(aux->pid, 15);
	}
	free(aux);
	return 0;
}
