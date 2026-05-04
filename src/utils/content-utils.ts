import I18nKey from '@i18n/i18nKey'
import { DEFAULT_LANG, LANGS, type Lang } from '@i18n/langs'
import { i18n } from '@i18n/translation'
import { getCollection } from 'astro:content'

export async function getSortedPosts(lang: Lang = DEFAULT_LANG) {
  const collectionName = LANGS[lang].postsCollection
  const allBlogPosts = await getCollection(collectionName, ({ data }) => {
    return import.meta.env.PROD ? data.draft !== true : true
  })
  const sorted = allBlogPosts.sort((a, b) => {
    const dateA = new Date(a.data.published)
    const dateB = new Date(b.data.published)
    return dateA > dateB ? -1 : 1
  })

  for (let i = 1; i < sorted.length; i++) {
    sorted[i].data.nextSlug = sorted[i - 1].slug
    sorted[i].data.nextTitle = sorted[i - 1].data.title
  }
  for (let i = 0; i < sorted.length - 1; i++) {
    sorted[i].data.prevSlug = sorted[i + 1].slug
    sorted[i].data.prevTitle = sorted[i + 1].data.title
  }

  return sorted
}

export type Tag = {
  name: string
  count: number
}

export async function getTagList(lang: Lang = DEFAULT_LANG): Promise<Tag[]> {
  const collectionName = LANGS[lang].postsCollection
  const allBlogPosts = await getCollection(collectionName, ({ data }) => {
    return import.meta.env.PROD ? data.draft !== true : true
  })

  const countMap: { [key: string]: number } = {}
  allBlogPosts.map(post => {
    post.data.tags.map((tag: string) => {
      if (!countMap[tag]) countMap[tag] = 0
      countMap[tag]++
    })
  })

  const keys: string[] = Object.keys(countMap).sort((a, b) => {
    return a.toLowerCase().localeCompare(b.toLowerCase())
  })

  return keys.map(key => ({ name: key, count: countMap[key] }))
}

export type Category = {
  name: string
  count: number
}

export async function getCategoryList(
  lang: Lang = DEFAULT_LANG,
): Promise<Category[]> {
  const collectionName = LANGS[lang].postsCollection
  const allBlogPosts = await getCollection(collectionName, ({ data }) => {
    return import.meta.env.PROD ? data.draft !== true : true
  })
  const count: { [key: string]: number } = {}
  allBlogPosts.map(post => {
    if (!post.data.category) {
      const ucKey = i18n(I18nKey.uncategorized, lang)
      count[ucKey] = count[ucKey] ? count[ucKey] + 1 : 1
      return
    }
    count[post.data.category] = count[post.data.category]
      ? count[post.data.category] + 1
      : 1
  })

  const lst = Object.keys(count).sort((a, b) => {
    return a.toLowerCase().localeCompare(b.toLowerCase())
  })

  const ret: Category[] = []
  for (const c of lst) {
    ret.push({ name: c, count: count[c] })
  }
  return ret
}
