import { LinkPreset, type NavBarLink } from '@/types/config'
import I18nKey from '@i18n/i18nKey'
import { DEFAULT_LANG, type Lang } from '@i18n/langs'
import { i18n } from '@i18n/translation'

export function getLinkPresets(
  lang: Lang = DEFAULT_LANG,
): { [key in LinkPreset]: NavBarLink } {
  return {
    [LinkPreset.Home]: {
      name: i18n(I18nKey.home, lang),
      url: '/',
    },
    [LinkPreset.About]: {
      name: i18n(I18nKey.about, lang),
      url: '/about/',
    },
    [LinkPreset.Archive]: {
      name: i18n(I18nKey.archive, lang),
      url: '/archive/',
    },
  }
}
