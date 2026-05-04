import i18nKey from '@i18n/i18nKey'
import { LANGS, type Lang, DEFAULT_LANG, getLangFromPath } from '@i18n/langs'
import { i18n } from '@i18n/translation'

export function pathsEqual(path1: string, path2: string) {
  const normalizedPath1 = path1.replace(/^\/|\/$/g, '').toLowerCase()
  const normalizedPath2 = path2.replace(/^\/|\/$/g, '').toLowerCase()
  return normalizedPath1 === normalizedPath2
}

function joinUrl(...parts: string[]): string {
  const joined = parts.join('/')
  return joined.replace(/\/+/g, '/')
}

function withLangPrefix(path: string, lang: Lang): string {
  const prefix = LANGS[lang].urlPrefix
  if (!prefix) return path
  const cleanPath = path.startsWith('/') ? path : `/${path}`
  return `/${prefix}${cleanPath}`
}

export function getPostUrlBySlug(
  slug: string,
  lang: Lang = DEFAULT_LANG,
): string {
  return url(withLangPrefix(`/posts/${slug}/`, lang))
}

export function getCategoryUrl(
  category: string,
  lang: Lang = DEFAULT_LANG,
): string {
  if (category === i18n(i18nKey.uncategorized, lang))
    return url(withLangPrefix('/archive/category/uncategorized/', lang))
  return url(withLangPrefix(`/archive/category/${category}/`, lang))
}

export function getDir(path: string): string {
  const lastSlashIndex = path.lastIndexOf('/')
  if (lastSlashIndex < 0) {
    return '/'
  }
  return path.substring(0, lastSlashIndex + 1)
}

export function url(path: string) {
  return joinUrl('', import.meta.env.BASE_URL, path)
}

export function localizedUrl(path: string, lang: Lang = DEFAULT_LANG): string {
  return url(withLangPrefix(path, lang))
}

export function getLangFromAstro(astroUrl: { pathname: string }): Lang {
  const base = import.meta.env.BASE_URL.replace(/\/+$/, '')
  let pathname = astroUrl.pathname
  if (base && pathname.startsWith(base)) {
    pathname = pathname.slice(base.length) || '/'
  }
  return getLangFromPath(pathname)
}
