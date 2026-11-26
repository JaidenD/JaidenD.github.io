/*
 * blog.js: turns the Markdown + LaTeX files in /posts into pages.
 *
 * You shouldn't need to edit this file to write posts. The two things you
 * might want to tweak live right at the top:
 *   MATH_MACROS:  shortcuts available in every post, e.g. $\R^n$.
 *   ENVIRONMENTS: the "::: theorem" / "::: proof" style blocks.
 */
(function () {
  'use strict';

  // ---------------------------------------------------------------------------
  // Settings
  // ---------------------------------------------------------------------------

  /** LaTeX macros available in every post. [body, number of arguments] for macros with arguments. */
  const MATH_MACROS = {
    R: '\\mathbb{R}',
    C: '\\mathbb{C}',
    Q: '\\mathbb{Q}',
    Z: '\\mathbb{Z}',
    N: '\\mathbb{N}',
    F: '\\mathbb{F}',
    eps: '\\varepsilon',
    abs: ['\\left\\lvert #1 \\right\\rvert', 1],
    norm: ['\\left\\lVert #1 \\right\\rVert', 1],
    inner: ['\\left\\langle #1, #2 \\right\\rangle', 2],
    Hom: '\\operatorname{Hom}',
    End: '\\operatorname{End}',
    Aut: '\\operatorname{Aut}',
    Gal: '\\operatorname{Gal}',
    tr: '\\operatorname{tr}',
    im: '\\operatorname{im}',
    id: '\\mathrm{id}',
  };

  /**
   * Theorem-like blocks. Write them as
   *
   *   ::: theorem Optional name
   *   Statement...
   *   :::
   *
   * Numbered environments share one counter, like amsthm's default.
   */
  const ENVIRONMENTS = {
    theorem: { label: 'Theorem', numbered: true, style: 'plain' },
    lemma: { label: 'Lemma', numbered: true, style: 'plain' },
    proposition: { label: 'Proposition', numbered: true, style: 'plain' },
    corollary: { label: 'Corollary', numbered: true, style: 'plain' },
    conjecture: { label: 'Conjecture', numbered: true, style: 'plain' },
    claim: { label: 'Claim', numbered: false, style: 'plain' },
    definition: { label: 'Definition', numbered: true, style: 'definition' },
    example: { label: 'Example', numbered: true, style: 'example' },
    exercise: { label: 'Exercise', numbered: true, style: 'example' },
    problem: { label: 'Problem', numbered: true, style: 'example' },
    axiom: { label: 'Axiom', numbered: false, style: 'definition' },
    principle: { label: 'Principle', numbered: false, style: 'definition' },
    remark: { label: 'Remark', numbered: false, style: 'remark' },
    note: { label: 'Note', numbered: false, style: 'remark' },
    proof: { label: 'Proof', numbered: false, style: 'proof' },
    solution: { label: 'Solution', numbered: false, style: 'proof' },
  };
  const ENV_ALIASES = {
    thm: 'theorem', lem: 'lemma', prop: 'proposition', cor: 'corollary',
    def: 'definition', defn: 'definition', ex: 'example', rem: 'remark',
  };

  const CDN = 'https://cdn.jsdelivr.net/npm/';
  const MATHJAX_URL = CDN + 'mathjax@3.2.2/es5/tex-chtml.js';
  const HLJS_URL = CDN + '@highlightjs/cdn-assets@11.11.1/highlight.min.js';

  const POSTS_DIR = 'posts/';
  const INDEX_URL = POSTS_DIR + 'index.json';
  const SITE_NAME = 'Jaiden Dahlke';
  const IS_LOCAL = ['localhost', '127.0.0.1', '[::1]', '::1'].includes(location.hostname);

  // ---------------------------------------------------------------------------
  // Small helpers
  // ---------------------------------------------------------------------------

  function escapeHtml(s) {
    return String(s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;')
      .replace(/"/g, '&quot;');
  }

  async function fetchText(url) {
    const res = await fetch(url, { cache: IS_LOCAL ? 'no-store' : 'no-cache' });
    if (!res.ok) throw new Error(res.status + ' ' + res.statusText);
    return res.text();
  }

  const loadedScripts = {};
  function loadScript(src) {
    if (!loadedScripts[src]) {
      loadedScripts[src] = new Promise((resolve, reject) => {
        const s = document.createElement('script');
        s.src = src;
        s.async = true;
        s.onload = resolve;
        s.onerror = () => reject(new Error('Could not load ' + src));
        document.head.appendChild(s);
      });
    }
    return loadedScripts[src];
  }

  function safeDecode(s) {
    try { return decodeURIComponent(s); } catch (e) { return s; }
  }

  function slugify(text) {
    return String(text)
      .toLowerCase()
      .normalize('NFKD')
      .replace(/[\u0300-\u036f]/g, '')
      .replace(/[^a-z0-9]+/g, '-')
      .replace(/^-+|-+$/g, '');
  }

  // ---------------------------------------------------------------------------
  // Front matter: a "---" block of "key: value" lines at the top of a post.
  // (tools/site.py parses exactly the same format.)
  // ---------------------------------------------------------------------------

  function parseScalar(raw) {
    const v = raw.trim();
    if ((v.startsWith('"') && v.endsWith('"') && v.length > 1) ||
        (v.startsWith("'") && v.endsWith("'") && v.length > 1)) return v.slice(1, -1);
    if (v.startsWith('[') && v.endsWith(']')) {
      return v.slice(1, -1).split(',').map((x) => parseScalar(x)).filter((x) => x !== '');
    }
    if (/^(true|yes)$/i.test(v)) return true;
    if (/^(false|no)$/i.test(v)) return false;
    return v;
  }

  function parseFrontMatter(src) {
    src = src.replace(/^\uFEFF/, '').replace(/\r\n?/g, '\n');
    const m = /^---[ \t]*\n([\s\S]*?)\n---[ \t]*(?:\n|$)/.exec(src);
    if (!m) return { meta: {}, body: src };
    const meta = {};
    let listKey = null;
    for (const line of m[1].split('\n')) {
      if (!line.trim() || line.trim().startsWith('#')) continue;
      const item = /^\s+-\s+(.*)$/.exec(line) || /^-\s+(.*)$/.exec(line);
      if (item && listKey) {
        meta[listKey].push(parseScalar(item[1]));
        continue;
      }
      const kv = /^([A-Za-z_][\w-]*)\s*:\s*(.*)$/.exec(line);
      if (!kv) continue;
      const key = kv[1].toLowerCase();
      if (kv[2].trim() === '') {
        meta[key] = [];
        listKey = key;
      } else {
        meta[key] = parseScalar(kv[2]);
        listKey = null;
      }
    }
    return { meta, body: src.slice(m[0].length) };
  }

  function asList(v) {
    if (Array.isArray(v)) return v.map(String).filter(Boolean);
    if (typeof v === 'string' && v.trim()) return v.split(',').map((s) => s.trim()).filter(Boolean);
    return [];
  }

  // ---------------------------------------------------------------------------
  // Math protection.
  //
  // Markdown and LaTeX disagree about characters like _ * \ and |. So before
  // the Markdown parser sees the text, every piece of math is swapped for a
  // placeholder. Afterwards the placeholders are swapped back and MathJax
  // typesets them. Recognised:
  //   $...$  \(...\)                    inline
  //   $$...$$  \[...\]  ```math fences  display
  //   \begin{align} ... \end{align}     (and equation, gather, multline, ...)
  //   \ref{...}, \eqref{...}            cross-references outside math
  // Code spans and code blocks are left alone, and \$ is a literal dollar sign.
  // ---------------------------------------------------------------------------

  const MATH_ENVS = [
    'equation', 'align', 'alignat', 'gather', 'multline', 'flalign', 'eqnarray', 'CD',
  ];
  const ENV_RE = new RegExp('\\\\begin\\{((?:' + MATH_ENVS.join('|') + ')\\*?)\\}', 'y');
  const REF_RE = /\\(?:eq)?ref\{[^{}\n]*\}/y;
  const FENCE_RE = /[ \t]*(`{3,}|~{3,})[ \t]*([^\n]*)/y;
  const BLANK_LINE_RE = /\n[ \t]*\n/g;
  const TOKEN_OPEN = '\uE000';
  const TOKEN_CLOSE = '\uE001';

  function paragraphEnd(src, from) {
    BLANK_LINE_RE.lastIndex = from;
    const m = BLANK_LINE_RE.exec(src);
    return m ? m.index : src.length;
  }

  function isSpace(c) {
    return c === undefined || /\s/.test(c);
  }

  // Closing $ of inline math. Like Pandoc: math stays within one paragraph,
  // "$x$" needs no space just inside the dollars, and the closing $ is not
  // followed by a digit. So "$5 and $10" and "$5 for $x$" stay as text.
  // "$ x $" (spaces on both sides) is accepted too.
  function findInlineClose(src, start) {
    if (start >= src.length || src[start] === '$') return -1;
    const padded = isSpace(src[start]);
    const limit = paragraphEnd(src, start);
    for (let j = start; j < limit; j++) {
      const c = src[j];
      if (c === '\\') { j++; continue; }
      if (c === '$') {
        const ok = isSpace(src[j - 1]) === padded && !/[0-9]/.test(src[j + 1] || '');
        return ok ? j : -1;
      }
    }
    return -1;
  }

  function findClose(src, close, from) {
    for (let j = from; j < src.length; j++) {
      if (src[j] === '\\' && close[0] !== '\\') { j++; continue; }
      if (src.startsWith(close, j)) return j;
    }
    return -1;
  }

  function protectMath(src) {
    const math = [];
    const out = [];
    let last = 0;
    let i = 0;
    const n = src.length;
    const emit = (end, tex, kind) => {
      out.push(src.slice(last, i));
      math.push({ tex, kind });
      out.push(TOKEN_OPEN + (math.length - 1) + TOKEN_CLOSE);
      i = last = end;
    };

    while (i < n) {
      const ch = src[i];

      // Fenced code block (``` or ~~~): skip it, unless it is a ```math block.
      if (i === 0 || src[i - 1] === '\n') {
        FENCE_RE.lastIndex = i;
        const f = FENCE_RE.exec(src);
        if (f && !(f[1][0] === '`' && f[2].includes('`'))) {
          const fence = f[1];
          const info = f[2].trim().toLowerCase();
          let lineStart = src.indexOf('\n', i);
          let closeStart = -1;
          let closeEnd = n;
          while (lineStart !== -1) {
            const lineEnd = src.indexOf('\n', lineStart + 1);
            const line = src.slice(lineStart + 1, lineEnd === -1 ? n : lineEnd);
            const c = /^[ \t]*(`{3,}|~{3,})[ \t]*$/.exec(line);
            if (c && c[1][0] === fence[0] && c[1].length >= fence.length) {
              closeStart = lineStart;
              closeEnd = lineEnd === -1 ? n : lineEnd;
              break;
            }
            lineStart = lineEnd;
          }
          if (info === 'math' && closeStart !== -1) {
            const bodyStart = src.indexOf('\n', i);
            emit(closeEnd, src.slice(bodyStart + 1, closeStart), 'display');
          } else {
            i = closeEnd;
          }
          continue;
        }
      }

      if (ch === '`') {
        // Inline code: a run of backticks closed by a run of the same length.
        let k = i;
        while (src[k] === '`') k++;
        const run = k - i;
        const limit = paragraphEnd(src, k);
        let j = k;
        let end = -1;
        while (j < limit) {
          const p = src.indexOf('`', j);
          if (p === -1 || p >= limit) break;
          let q = p;
          while (src[q] === '`') q++;
          if (q - p === run) { end = q; break; }
          j = q;
        }
        i = end === -1 ? k : end;
        continue;
      }

      if (ch === '\\') {
        const next = src[i + 1];
        if (next === '[' || next === '(') {
          const close = next === '[' ? '\\]' : '\\)';
          const end = findClose(src, close, i + 2);
          if (end !== -1) {
            emit(end + 2, src.slice(i + 2, end), next === '[' ? 'display' : 'inline');
            continue;
          }
        } else if (next === 'b') {
          ENV_RE.lastIndex = i;
          const m = ENV_RE.exec(src);
          if (m) {
            const endTag = '\\end{' + m[1] + '}';
            const end = src.indexOf(endTag, ENV_RE.lastIndex);
            if (end !== -1) {
              emit(end + endTag.length, src.slice(i, end + endTag.length), 'env');
              continue;
            }
          }
        } else if (next === 'r' || next === 'e') {
          REF_RE.lastIndex = i;
          const m = REF_RE.exec(src);
          if (m) {
            emit(REF_RE.lastIndex, m[0], 'inline');
            continue;
          }
        }
        i += 2; // an escaped character such as \$ or \*; Markdown handles it
        continue;
      }

      if (ch === '$') {
        if (src[i + 1] === '$') {
          const end = findClose(src, '$$', i + 2);
          if (end !== -1) {
            emit(end + 2, src.slice(i + 2, end), 'display');
            continue;
          }
          i += 2;
          continue;
        }
        const end = findInlineClose(src, i + 1);
        if (end !== -1) {
          emit(end + 1, src.slice(i + 1, end), 'inline');
          continue;
        }
      }
      i++;
    }
    out.push(src.slice(last));
    return { text: out.join(''), math };
  }

  // Math inside an HTML attribute (image alt text, a title) can't become markup,
  // so it goes back as plain TeX with its dollar signs.
  function mathAsText(m) {
    if (!m) return '';
    if (m.kind === 'env') return m.tex;
    return m.kind === 'inline' ? '$' + m.tex + '$' : '$$' + m.tex + '$$';
  }

  const DEFINITIONS_RE = /^\s*\\(?:(?:re)?newcommand|def|let|DeclareMathOperator|newenvironment)\b/;
  const TOKEN_RE = new RegExp(TOKEN_OPEN + '(\\d+)' + TOKEN_CLOSE, 'g');
  const PARA_TOKEN_RE = new RegExp('<p>\\s*' + TOKEN_OPEN + '(\\d+)' + TOKEN_CLOSE + '\\s*</p>', 'g');

  function mathHtml(m, alone) {
    if (!m) return '';
    if (m.kind === 'inline') {
      const span = '<span class="math tex2jax_process">' + escapeHtml('\\(' + m.tex + '\\)') + '</span>';
      if (!alone) return span;
      // A paragraph holding nothing but \newcommand-style definitions: keep it invisible.
      if (DEFINITIONS_RE.test(m.tex)) return '<div class="math-defs tex2jax_process">' + escapeHtml('\\(' + m.tex + '\\)') + '</div>';
      return '<p>' + span + '</p>';
    }
    const body = m.kind === 'env' ? m.tex : '\\[' + m.tex + '\\]';
    const cls = alone && DEFINITIONS_RE.test(m.tex) ? 'math-defs' : 'math-display';
    const tag = alone ? 'div' : 'span';
    return '<' + tag + ' class="' + cls + ' tex2jax_process">' + escapeHtml(body) + '</' + tag + '>';
  }

  function restoreMath(html, math) {
    const inAttr = (attr) => attr.replace(TOKEN_RE, (_, k) => escapeHtml(mathAsText(math[+k])));
    return html
      .replace(/="[^"]*"/g, inAttr)
      .replace(/='[^']*'/g, inAttr)
      .replace(PARA_TOKEN_RE, (_, k) => mathHtml(math[+k], true))
      .replace(TOKEN_RE, (_, k) => mathHtml(math[+k], false));
  }

  /** Plain text that may contain $math$ (titles, summaries) -> safe HTML. */
  function textWithMath(str) {
    const { text, math } = protectMath(String(str || ''));
    return {
      html: restoreMath(escapeHtml(text.replace(/\\\$/g, '$')), math),
      hasMath: math.length > 0,
    };
  }

  // ---------------------------------------------------------------------------
  // MathJax (loaded only when a page actually contains math)
  // ---------------------------------------------------------------------------

  let mathJaxPromise = null;

  function loadMathJax() {
    if (mathJaxPromise) return mathJaxPromise;
    mathJaxPromise = new Promise((resolve, reject) => {
      window.MathJax = {
        loader: { load: ['[tex]/mathtools'] },
        tex: {
          inlineMath: [['\\(', '\\)']],
          displayMath: [['\\[', '\\]']],
          processEscapes: false,
          processEnvironments: true,
          processRefs: false,
          tags: 'ams',
          macros: MATH_MACROS,
          packages: { '[+]': ['mathtools'] },
        },
        options: {
          ignoreHtmlClass: 'tex2jax_ignore',
          processHtmlClass: 'tex2jax_process',
        },
        startup: {
          typeset: false,
          ready() {
            window.MathJax.startup.defaultReady();
            window.MathJax.startup.promise.then(() => resolve(window.MathJax));
          },
        },
      };
      loadScript(MATHJAX_URL).catch(reject);
    });
    return mathJaxPromise;
  }

  async function typeset(root) {
    const MJ = await loadMathJax();
    // MathJax sizes math to match the text font, so let the web fonts arrive first.
    if (document.fonts && document.fonts.ready) {
      await Promise.race([document.fonts.ready, new Promise((r) => setTimeout(r, 1500))]);
    }
    MJ.typesetClear();
    MJ.texReset();
    await MJ.typesetPromise([root]);
  }

  // ---------------------------------------------------------------------------
  // Markdown
  // ---------------------------------------------------------------------------

  function resolveEnv(info) {
    const [first, ...rest] = info.trim().split(/\s+/);
    const word = (first || '').toLowerCase();
    const key = ENVIRONMENTS[word] ? word : ENV_ALIASES[word];
    if (!key) return null;
    return { key, env: ENVIRONMENTS[key], title: rest.join(' ') };
  }

  function createMarkdown() {
    const md = window.markdownit({ html: true, linkify: false, typographer: true });
    // Keep smart quotes, but not the (c) -> © style replacements: "(c)" is a
    // perfectly good label for part (c) of an exercise.
    md.disable('replacements');
    if (window.markdownitFootnote) {
      md.use(window.markdownitFootnote);
      md.renderer.rules.footnote_caption = (tokens, idx) => {
        const n = Number(tokens[idx].meta.id + 1).toString();
        return tokens[idx].meta.subId > 0 ? n + ':' + tokens[idx].meta.subId : n;
      };
    }

    // Theorem numbers; reset before each render in renderMarkdown().
    const state = md.envState = { counter: 0 };

    if (window.markdownitContainer) {
      md.use(window.markdownitContainer, 'env', {
        validate: (params) => resolveEnv(params) !== null,
        render(tokens, idx) {
          const tok = tokens[idx];
          if (tok.nesting !== 1) return '</div>\n';
          const { key, env, title } = resolveEnv(tok.info);
          let head;
          if (env.style === 'proof') {
            head = escapeHtml(env.label) + (title ? ' ' + md.renderInline(title) : '');
          } else {
            head = escapeHtml(env.label) + (env.numbered ? ' ' + ++state.counter : '');
            if (title) head += ' <span class="env-title">(' + md.renderInline(title) + ')</span>';
          }
          return '<div class="env env-' + key + ' env-style-' + env.style + '">' +
            '<span class="env-head">' + head + '.</span>\n';
        },
      });
    }
    return md;
  }

  let markdown = null;
  function renderMarkdown(body) {
    if (!markdown) markdown = createMarkdown();
    const { text, math } = protectMath(body);
    markdown.envState.counter = 0;
    return restoreMath(markdown.render(text), math);
  }

  // Tidy-ups that are easier on the DOM than in the Markdown renderer.
  function enhance(root) {
    // Theorem heads run into the first paragraph, like in LaTeX.
    root.querySelectorAll('.env > .env-head').forEach((head) => {
      const next = head.nextElementSibling;
      if (next && next.tagName === 'P') {
        next.insertBefore(document.createTextNode(' '), next.firstChild);
        next.insertBefore(head, next.firstChild);
      } else {
        const p = document.createElement('p');
        head.replaceWith(p);
        p.appendChild(head);
      }
    });

    // End of proof.
    root.querySelectorAll('.env-style-proof').forEach((proof) => {
      const qed = document.createElement('span');
      qed.className = 'qed';
      qed.setAttribute('aria-label', 'end of proof');
      qed.textContent = '\u220E';
      let last = proof.lastElementChild;
      if (!last || last.tagName !== 'P') {
        last = document.createElement('p');
        last.className = 'qed-line';
        proof.appendChild(last);
      }
      last.appendChild(qed);
    });

    // Heading ids + anchor links.
    const used = new Set();
    root.querySelectorAll('h2, h3, h4').forEach((h, n) => {
      if (!h.id) {
        const clone = h.cloneNode(true);
        clone.querySelectorAll('.tex2jax_process').forEach((m) => m.remove());
        let id = slugify(clone.textContent) || 'section-' + (n + 1);
        while (used.has(id) || document.getElementById(id)) id += '-';
        h.id = id;
      }
      used.add(h.id);
      const a = document.createElement('a');
      a.className = 'heading-anchor';
      a.href = '#' + h.id;
      a.setAttribute('aria-label', 'Link to this section');
      a.textContent = '#';
      h.appendChild(a);
    });

    // Lone images with a title become figures with a caption: ![alt](src "Caption").
    root.querySelectorAll('p > img:only-child').forEach((img) => {
      const p = img.parentElement;
      if (p.textContent.trim()) return;
      const fig = document.createElement('figure');
      fig.appendChild(img);
      if (img.title) {
        const cap = document.createElement('figcaption');
        cap.innerHTML = textWithMath(img.title).html;
        fig.appendChild(cap);
        img.removeAttribute('title');
      }
      p.replaceWith(fig);
    });

    // Paths in a post are relative to the .md file, so "images/plot.png" means
    // posts/images/plot.png, and a link to "other-post.md" opens that post.
    root.querySelectorAll('img[src], a[href], video[src], source[src], script[src]').forEach((el) => {
      const attr = el.hasAttribute('href') ? 'href' : 'src';
      const url = el.getAttribute(attr);
      if (!url || /^(?:[a-z][a-z0-9+.-]*:|\/|#|\?)/i.test(url)) return;
      const post = /^(?:\.\/)?([^/?#]+)\.md(#.*)?$/.exec(url);
      el.setAttribute(attr, post
        ? 'post.html?p=' + encodeURIComponent(safeDecode(post[1])) + (post[2] || '')
        : POSTS_DIR + url);
    });

    // Scripts in a post (say, for an interactive figure) only run if re-created.
    root.querySelectorAll('script').forEach((old) => {
      const script = document.createElement('script');
      for (const a of old.attributes) script.setAttribute(a.name, a.value);
      script.async = false;
      script.textContent = old.textContent;
      old.replaceWith(script);
    });

    // Wide tables scroll sideways instead of breaking the layout on phones.
    root.querySelectorAll('table').forEach((t) => {
      const wrap = document.createElement('div');
      wrap.className = 'table-wrap';
      t.replaceWith(wrap);
      wrap.appendChild(t);
    });
  }

  function buildToc(root, nav) {
    const heads = root.querySelectorAll('h2, h3');
    if (heads.length < 2) { nav.hidden = true; return; }
    let html = '<p class="toc-title">Contents</p><ol>';
    let open = false;
    heads.forEach((h) => {
      const clone = h.cloneNode(true);
      clone.querySelectorAll('.heading-anchor').forEach((a) => a.remove());
      const link = '<a href="#' + h.id + '">' + clone.innerHTML + '</a>';
      if (h.tagName === 'H3') {
        if (!open) { html += '<ol>'; open = true; }
        html += '<li>' + link + '</li>';
      } else {
        if (open) { html += '</ol>'; open = false; }
        html += '<li>' + link;
      }
    });
    if (open) html += '</ol>';
    nav.innerHTML = html + '</ol>';
    nav.hidden = false;
  }

  async function highlightCode(root) {
    const blocks = root.querySelectorAll('pre code[class*="language-"]');
    if (!blocks.length) return;
    try {
      await loadScript(HLJS_URL);
      window.hljs.configure({ ignoreUnescapedHTML: true });
      blocks.forEach((b) => {
        const lang = (/language-(\S+)/.exec(b.className) || [])[1];
        if (lang && window.hljs.getLanguage(lang)) window.hljs.highlightElement(b);
      });
    } catch (e) { /* plain code blocks are fine too */ }
  }

  // ---------------------------------------------------------------------------
  // Pages
  // ---------------------------------------------------------------------------

  const MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
    'August', 'September', 'October', 'November', 'December'];

  /** "2026-06-18" -> "18 June 2026" (or "June 2026" with short: true); "2026-06" -> "June 2026". */
  function formatDate(value, short) {
    const m = /^(\d{4})-(\d{1,2})(?:-(\d{1,2}))?$/.exec(String(value || '').trim());
    if (!m || +m[2] < 1 || +m[2] > 12) return String(value || '');
    const month = MONTHS[+m[2] - 1] + ' ' + m[1];
    return m[3] && !short ? +m[3] + ' ' + month : month;
  }

  function formatTags(tags) {
    return asList(tags).map((t) => '<span class="tag">#' + escapeHtml(t) + '</span>').join(' ');
  }

  function plainTitle(title) {
    return String(title).replace(/\$+([^$]*)\$+/g, '$1');
  }

  async function renderPost(src, slug, els) {
    const { meta, body } = parseFrontMatter(src);
    const title = meta.title || slug;
    const date = meta.date ? String(meta.date) : '';

    const titleHtml = textWithMath(title);
    els.title.innerHTML = titleHtml.html;
    document.title = plainTitle(title) + ' \u2014 ' + SITE_NAME;

    const metaBits = [];
    if (meta.author) metaBits.push('<span class="author">by ' + escapeHtml(meta.author) + '</span>');
    if (date) metaBits.push('<time datetime="' + escapeHtml(date) + '">' + escapeHtml(formatDate(date)) + '</time>');
    if (meta.status) metaBits.push('<span class="status-badge">' + escapeHtml(meta.status) + '</span>');
    if (meta.draft === true) metaBits.push('<span class="draft-badge">draft</span>');
    const tags = formatTags(meta.tags);
    if (tags) metaBits.push(tags);
    els.meta.innerHTML = metaBits.join('<span class="sep"> \u00B7 </span>');

    try {
      els.body.innerHTML = renderMarkdown(body);
    } catch (e) {
      // The Markdown library didn't load (offline, CDN down): still show the text.
      console.error(e);
      els.body.innerHTML = '<pre class="raw-post">' + escapeHtml(body) + '</pre>';
      return;
    }
    enhance(els.body);
    if (meta.toc === true) buildToc(els.body, els.toc);
    else els.toc.hidden = true;

    const jobs = [highlightCode(els.body)];
    if (els.article.querySelector('.tex2jax_process')) jobs.push(typeset(els.article));
    await Promise.all(jobs);
  }

  function reveal(article, loading) {
    article.classList.remove('is-loading');
    if (loading) loading.hidden = true;
  }

  async function renderPostPage() {
    const els = {
      article: document.getElementById('post'),
      title: document.getElementById('post-title'),
      meta: document.getElementById('post-meta'),
      toc: document.getElementById('post-toc'),
      body: document.getElementById('post-body'),
    };
    const loading = document.getElementById('loading');
    const slug = new URLSearchParams(location.search).get('p') || '';
    const safe = slug && !/[\\/]|^\./.test(slug);
    const url = POSTS_DIR + encodeURIComponent(slug) + '.md';

    let src;
    try {
      if (!safe) throw new Error('no post given');
      src = await fetchText(url);
    } catch (e) {
      els.title.textContent = 'Post not found';
      els.meta.textContent = '';
      els.body.innerHTML = '<p>There is no post called <code>' + escapeHtml(slug || '(nothing)') +
        '</code>. <a href="./">Back to the home page.</a></p>' +
        (location.protocol === 'file:' ? fileProtocolHelp() : '');
      document.title = 'Not found \u2014 ' + SITE_NAME;
      reveal(els.article, loading);
      return;
    }

    // Show the page even if MathJax is slow or blocked.
    const failsafe = setTimeout(() => reveal(els.article, loading), 4000);
    try {
      await renderPost(src, slug, els);
    } catch (e) {
      console.error(e);
    }
    clearTimeout(failsafe);
    reveal(els.article, loading);
    if (location.hash) {
      const target = document.getElementById(decodeURIComponent(location.hash.slice(1)));
      if (target) target.scrollIntoView();
    }

    // Live preview while writing: re-render whenever the file changes.
    if (IS_LOCAL) {
      let lastSrc = src;
      let busy = false;
      setInterval(async () => {
        if (busy || document.hidden) return;
        busy = true;
        try {
          const next = await fetchText(url);
          if (next !== lastSrc) {
            lastSrc = next;
            const y = window.scrollY;
            await renderPost(next, slug, els);
            window.scrollTo(0, y);
          }
        } catch (e) { /* file mid-save or server stopped */ }
        busy = false;
      }, 1000);
    }
  }

  function fileProtocolHelp() {
    return '<p class="help">This page was opened straight from disk, so the browser won\u2019t let it read the ' +
      'post files. Preview the site with <code>python3 tools/site.py serve</code> and open ' +
      '<a href="http://localhost:8000/">http://localhost:8000/</a> instead.</p>';
  }

  async function renderIndexPage() {
    const list = document.getElementById('post-list');
    let posts;
    try {
      posts = JSON.parse(await fetchText(INDEX_URL));
    } catch (e) {
      list.innerHTML = '<li class="empty">Couldn\u2019t load the list of posts.' +
        (location.protocol === 'file:' ? fileProtocolHelp() : '') + '</li>';
      return;
    }
    posts = (Array.isArray(posts) ? posts : [])
      .filter((p) => p && p.slug)
      .sort((a, b) => String(b.date || '').localeCompare(String(a.date || '')));
    if (!posts.length) {
      list.innerHTML = '<li class="empty">Nothing here yet.</li>';
      return;
    }
    let needsMath = false;
    list.innerHTML = posts.map((p) => {
      const title = textWithMath(p.title || p.slug);
      const summary = p.summary ? textWithMath(p.summary) : null;
      needsMath = needsMath || title.hasMath || (summary && summary.hasMath);
      const bits = [];
      if (p.status) bits.push('<span class="status-badge">' + escapeHtml(p.status) + '</span>');
      if (p.draft) bits.push('<span class="draft-badge">draft</span>');
      if (p.date) bits.push('<time datetime="' + escapeHtml(p.date) + '">' + escapeHtml(formatDate(p.date, true)) + '</time>');
      return '<li class="post-item' + (p.draft ? ' is-draft' : '') + '">' +
        '<a class="post-item-title" href="post.html?p=' + encodeURIComponent(p.slug) + '">' + title.html + '</a>' +
        (bits.length ? ' <span class="post-item-meta">' + bits.join(' ') + '</span>' : '') +
        (summary ? '<p class="post-item-summary">' + summary.html + '</p>' : '') +
        '</li>';
    }).join('');
    if (needsMath) typeset(list).catch((e) => console.error(e));
  }

  function main() {
    const page = document.body.dataset.page;
    if (page === 'post') renderPostPage();
    else if (page === 'home') renderIndexPage();
  }

  // Exposed for tools/tests; the site itself doesn't need these.
  window.Blog = { parseFrontMatter, protectMath, restoreMath, renderMarkdown };

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', main);
  else main();
})();
