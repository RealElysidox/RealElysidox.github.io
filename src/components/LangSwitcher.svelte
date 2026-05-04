<script lang="ts">
import Icon from '@iconify/svelte'

type Lang = 'en' | 'zh_CN' | 'ja'

interface LangInfo {
  code: Lang
  urlPrefix: string
  label: string
  longLabel: string
}

const LANG_LIST: LangInfo[] = [
  { code: 'en', urlPrefix: '', label: 'EN', longLabel: 'English' },
  { code: 'zh_CN', urlPrefix: 'zh-cn', label: '中', longLabel: '简体中文' },
  { code: 'ja', urlPrefix: 'ja', label: '日', longLabel: '日本語' },
]

const BASE_URL = (import.meta.env.BASE_URL || '/').replace(/\/+$/, '')

export let lang: Lang = 'en'

const current: LangInfo = LANG_LIST.find(l => l.code === lang) ?? LANG_LIST[0]

function stripBase(pathname: string): string {
  if (BASE_URL && pathname.startsWith(BASE_URL)) {
    return pathname.slice(BASE_URL.length) || '/'
  }
  return pathname
}

function stripLang(pathname: string): string {
  for (const l of LANG_LIST) {
    if (!l.urlPrefix) continue
    const prefix = `/${l.urlPrefix}`
    if (pathname === prefix) return '/'
    if (pathname.startsWith(`${prefix}/`)) return pathname.slice(prefix.length)
  }
  return pathname
}

function withBase(path: string): string {
  if (!BASE_URL) return path
  return `${BASE_URL}${path}`.replace(/\/+/g, '/') || '/'
}

function buildTarget(target: LangInfo): string {
  if (typeof window === 'undefined') return '#'
  const stripped = stripLang(stripBase(window.location.pathname))
  // For posts, slugs differ between languages — fall back to homepage
  const isPost = /^\/posts(\/|$)/.test(stripped)
  const finalPath = isPost
    ? '/'
    : stripped.startsWith('/')
      ? stripped
      : `/${stripped}`
  const prefixed = target.urlPrefix
    ? `/${target.urlPrefix}${finalPath === '/' ? '/' : finalPath}`
    : finalPath
  return withBase(prefixed)
}

function switchTo(target: LangInfo) {
  if (target.code === lang) {
    hidePanel()
    return
  }
  window.location.href = buildTarget(target)
}

function showPanel() {
  const panel = document.querySelector('#lang-panel')
  panel?.classList.remove('float-panel-closed')
}

function hidePanel() {
  const panel = document.querySelector('#lang-panel')
  panel?.classList.add('float-panel-closed')
}
</script>

<div class="relative z-50" role="menu" tabindex="-1" on:mouseleave={hidePanel}>
    <button
        aria-label="Language"
        role="menuitem"
        class="btn-plain scale-animation rounded-lg h-11 w-11 active:scale-90 flex items-center justify-center"
        on:click={showPanel}
        on:mouseenter={showPanel}
    >
        <Icon icon="material-symbols:translate" class="text-[1.25rem]"></Icon>
    </button>

    <div id="lang-panel" class="hidden lg:block absolute transition float-panel-closed top-11 -right-2 pt-5">
        <div class="card-base float-panel p-2">
            {#each LANG_LIST as l}
                <button
                    class="flex transition whitespace-nowrap items-center justify-start w-full btn-plain scale-animation rounded-lg h-9 px-3 font-medium active:scale-95 mb-0.5 last:mb-0"
                    class:current-lang-btn={l.code === current.code}
                    on:click={() => switchTo(l)}
                >
                    <span class="inline-block w-6 mr-2 text-center">{l.label}</span>
                    {l.longLabel}
                </button>
            {/each}
        </div>
    </div>
</div>

<style lang="css">
.current-lang-btn {
    background: var(--btn-plain-bg-hover);
    color: var(--primary);
}
</style>
