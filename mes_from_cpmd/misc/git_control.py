
def get_git_version():
        from commit_stuff.version_hash import commit_hash, commit_date, commit_message
        print("# Hello. I am from commit {}".format(commit_hash))
        print("# Commit Date: {}".format(commit_date))
        print("# Commit Message: {}".format(commit_message))

