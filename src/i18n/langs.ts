export type Lang = 'en' | 'zh_CN' | 'ja'

export interface LangInfo {
  code: Lang
  urlPrefix: string
  htmlLang: string
  label: string
  longLabel: string
  postsCollection: 'posts' | 'posts-zh-cn' | 'posts-ja'
  aboutEntry: 'about' | 'about-zh-cn' | 'about-ja'
}

export const LANGS: Record<Lang, LangInfo> = {
  en: {
    code: 'en',
    urlPrefix: '',
    htmlLang: 'en',
    label: 'EN',
    longLabel: 'English',
    postsCollection: 'posts',
    aboutEntry: 'about',
  },
  zh_CN: {
    code: 'zh_CN',
    urlPrefix: 'zh-cn',
    htmlLang: 'zh-CN',
    label: '中',
    longLabel: '简体中文',
    postsCollection: 'posts-zh-cn',
    aboutEntry: 'about-zh-cn',
  },
  ja: {
    code: 'ja',
    urlPrefix: 'ja',
    htmlLang: 'ja',
    label: '日',
    longLabel: '日本語',
    postsCollection: 'posts-ja',
    aboutEntry: 'about-ja',
  },
}

export const ALL_LANGS: Lang[] = ['en', 'zh_CN', 'ja']

export const DEFAULT_LANG: Lang = 'en'

export function getLangFromPath(pathname: string): Lang {
  const segments = pathname.replace(/^\/+|\/+$/g, '').split('/')
  const first = segments[0]?.toLowerCase()
  if (first === 'zh-cn') return 'zh_CN'
  if (first === 'ja') return 'ja'
  return 'en'
}

export function stripLangFromPath(pathname: string): string {
  const lang = getLangFromPath(pathname)
  if (lang === 'en') return pathname
  const prefix = `/${LANGS[lang].urlPrefix}`
  if (pathname === prefix) return '/'
  if (pathname.startsWith(`${prefix}/`)) return pathname.slice(prefix.length)
  return pathname
}
