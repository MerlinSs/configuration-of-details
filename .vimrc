set encoding=utf-8		" 打开文件时编码格式
set fileencoding=utf-8          " 在保存文件时，指定编码
set termencoding=utf-8          " 终端环境告诉vim使用编码
"set mouse=a                     " 启用鼠标
set shiftwidth=4                " 缩进的空格数
" plug issues.
call plug#begin()

Plug 'lervag/vimtex'
Plug 'SirVer/ultisnips'
Plug 'honza/vim-snippets'
Plug 'chusiang/vim-sdcv'

call plug#end()

" This is necessary for VimTeX to load properly. The "indent" is optional.
" Note that most plugin managers will do this automatically.
filetype plugin indent on

" This enables Vim's and neovim's syntax-related features. Without this, some
" VimTeX features will not work (see ":help vimtex-requirements" for more
" info).
syntax enable
let g:tex_flavor='latex'
let g:vimtex_view_method='general'
let g:vimtex_quickfix_mode=0
set conceallevel=1
let g:tex_conceal='abdmg'
" 对中文的支持
let g:Tex_CompileRule_pdf = 'xelatex -synctex=1 --interaction=nonstopmode $*'
let g:vimtex_compiler_latexmk_engines = {'_':'-xelatex'}
let g:vimtex_compiler_latexrun_engines ={'_':'xelatex'}

set conceallevel=2
let g:tex_conceal='abdmg'

"set number
"highlight LineNr cterm=bold ctermfg=7
"highlight LineNr ctermbg=237

nmap <leader>w :call SearchWord()<CR>		"查字典
