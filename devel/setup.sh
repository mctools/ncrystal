if [ "x${BASH_VERSION:-}" != "x" ]; then
    export PATH="$( cd -P "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/bin:${PATH}"
    type complete compgen 1>/dev/null 2>/dev/null
    if [ $? == 0 ]; then
        _ncdevtool_completions ()
        {
            COMPREPLY=($(compgen -W "$(ncdevtool --show-completion-list)" "${COMP_WORDS[1]}"))
        }
        complete -F _ncdevtool_completions ncdevtool
    fi
elif [ "x${ZSH_VERSION:-}" != "x" ]; then
    export PATH="${0:a:h}/bin:${PATH}"
else
    echo "ERROR: This script only supports bash or zsh for now. You will"
    echo "       have to manually add <reporoot>/devel/bin to your PATH"
fi
