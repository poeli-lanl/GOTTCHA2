#!/usr/bin/env python3

import argparse
import base64
import bisect
import gzip
import os
import json
import math
import re
import sys
import tempfile
from pathlib import Path
import pandas as pd
from collections import defaultdict

try:
    import minify_html
except ImportError:
    minify_html = None

HTML_ASSET_DIR = os.path.realpath(
    os.path.join(os.path.dirname(__file__), '..', 'data', 'html')
)

RESOURCE_SPECS = (
    {
        'marker': '<script src="/publicdata/js/d3.v7.min.js"></script>',
        'path': 'js/d3.v7.min.js',
        'kind': 'script',
    },
    {
        'marker': '<script src="/publicdata/js/vue.global.prod.js"></script>',
        'path': 'js/vue.global.prod.js',
        'kind': 'script',
    },
    {
        'marker': '<script src="/publicdata/js/primevue.min.js"></script>',
        'path': 'js/primevue.min.js',
        'kind': 'script',
    },
    {
        'marker': '<script src="/publicdata/js/aura.js"></script>',
        'path': 'js/aura.js',
        'kind': 'script',
    },
    {
        'marker': '<link href="/publicdata/css/bootstrap.min.css" rel="stylesheet">',
        'path': 'css/bootstrap.min.css',
        'kind': 'style',
    },
    {
        'marker': '<link rel="stylesheet" href="/publicdata/css/primeicons.css">',
        'path': 'css/primeicons.css',
        'kind': 'primeicons',
    },
    {
        'marker': '<script src="/publicdata/js/bootstrap.bundle.min.js"></script>',
        'path': 'js/bootstrap.bundle.min.js',
        'kind': 'script',
    },
)


def _read_resource(source):
    """Read a browser dependency from the bundled HTML asset directory."""
    try:
        with open(source, 'rb') as resource_file:
            return resource_file.read()
    except OSError as exc:
        raise RuntimeError(f"Could not read browser resource {source}: {exc}") from exc


def _embed_primeicons_font(css, stylesheet_source):
    """Replace PrimeIcons' multi-format font references with one embedded WOFF2 font."""
    font_source = os.path.join(
        os.path.dirname(stylesheet_source), 'fonts', 'primeicons.woff2'
    )
    font_data = base64.b64encode(_read_resource(font_source)).decode('ascii')
    embedded_font_face = (
        "@font-face{font-family:'primeicons';font-display:block;"
        f"src:url('data:font/woff2;base64,{font_data}') format('woff2');"
        "font-weight:normal;font-style:normal}"
    )
    css, substitutions = re.subn(
        r'@font-face\s*\{.*?\}',
        embedded_font_face,
        css,
        count=1,
        flags=re.DOTALL,
    )
    if substitutions != 1:
        raise ValueError(
            'PrimeIcons stylesheet does not contain the expected @font-face rule'
        )
    return css


def _embed_browser_resources(html_template):
    """Embed local browser dependencies as one deterministic compressed payload."""
    resources = []
    for spec in RESOURCE_SPECS:
        source = os.path.join(HTML_ASSET_DIR, spec['path'])
        content = _read_resource(source).decode('utf-8')
        content = re.sub(
            r'/\*[#@]\s*sourceMappingURL=.*?\*/', '', content, flags=re.DOTALL
        )
        content = re.sub(
            r'^\s*//[#@]\s*sourceMappingURL=.*$', '', content, flags=re.MULTILINE
        )
        if spec['kind'] == 'primeicons':
            content = _embed_primeicons_font(content, source)
        if spec['marker'] not in html_template:
            raise ValueError(f"HTML resource marker not found: {spec['marker']}")
        resources.append({
            'kind': 'script' if spec['kind'] == 'script' else 'style',
            'content': content,
        })

    resource_json = json.dumps(
        {'resources': resources},
        ensure_ascii=False,
        separators=(',', ':'),
    ).encode('utf-8')
    resource_payload = base64.b64encode(
        gzip.compress(resource_json, compresslevel=9, mtime=0)
    ).decode('ascii')
    style_count = sum(resource['kind'] == 'style' for resource in resources)

    loader = f"""<script>
window.browserResourcesReady = (async function() {{
  if (typeof DecompressionStream === 'undefined') {{
    throw new Error('This report requires a browser with DecompressionStream support.');
  }}

  const styleSlots = Array.from({{ length: {style_count} }}, function() {{
    const style = document.createElement('style');
    document.head.appendChild(style);
    return style;
  }});
  const binary = atob('{resource_payload}');
  const bytes = new Uint8Array(binary.length);
  for (let index = 0; index < binary.length; index += 1) {{
    bytes[index] = binary.charCodeAt(index);
  }}
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
  const payload = JSON.parse(await new Response(stream).text());
  let styleIndex = 0;

  for (const resource of payload.resources) {{
    if (resource.kind === 'style') {{
      styleSlots[styleIndex].textContent = resource.content;
      styleIndex += 1;
    }} else {{
      const script = document.createElement('script');
      script.textContent = resource.content;
      document.head.appendChild(script);
    }}
  }}
}})();
</script>"""

    html_content = html_template.replace(RESOURCE_SPECS[0]['marker'], loader, 1)
    for spec in RESOURCE_SPECS[1:]:
        html_content = html_content.replace(spec['marker'], '', 1)
    return html_content

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GOTTCHA2 Signature Coverage</title>
    <script src="/publicdata/js/d3.v7.min.js"></script>
    <script src="/publicdata/js/vue.global.prod.js"></script>
    <script src="/publicdata/js/primevue.min.js"></script>
    <script src="/publicdata/js/aura.js"></script>
    <link href="/publicdata/css/bootstrap.min.css" rel="stylesheet">
    <link rel="stylesheet" href="/publicdata/css/primeicons.css">
    <style>
        :root {
            --bg: #f6f8fb;
            --panel: #ffffff;
            --panel-soft: #f8fafc;
            --text: #172033;
            --muted: #697386;
            --line: #dbe3ef;
            --line-strong: #c2ccda;
            --accent: #2563eb;
            --accent-dark: #1d4ed8;
            --accent-soft: #dbeafe;
            --teal: #0f766e;
            --rose: #be3455;
            --shadow: 0 18px 45px rgba(31, 41, 55, 0.08);
        }

        * {
            box-sizing: border-box;
        }

        body {
            font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
            margin: 0;
            min-height: 100vh;
            background: var(--bg);
            color: var(--text);
        }

        button,
        select {
            font: inherit;
        }

        .app-shell {
            /* max-width: 1280px; */
            margin: 0 auto;
            padding: 28px;
            visibility: hidden;
        }

        #page-loader {
            align-items: center;
            background: var(--bg);
            display: flex;
            flex-direction: column;
            gap: 14px;
            inset: 0;
            justify-content: center;
            position: fixed;
            z-index: 2000;
        }

        .page-loader-spinner {
            animation: page-loader-spin 0.8s linear infinite;
            border: 4px solid var(--line);
            border-radius: 50%;
            border-top-color: var(--accent);
            height: 42px;
            width: 42px;
        }

        .page-loader-text {
            color: var(--muted);
            font-size: 13px;
            font-weight: 650;
        }

        @keyframes page-loader-spin {
            to {
                transform: rotate(360deg);
            }
        }

        .app-header {
            display: flex;
            align-items: flex-end;
            justify-content: space-between;
            gap: 24px;
            margin-bottom: 20px;
        }

        .eyebrow {
            color: var(--accent);
            font-size: 12px;
            font-weight: 750;
            margin: 0 0 6px;
            text-transform: uppercase;
        }

        h1 {
            color: var(--text);
            font-size: 30px;
            font-weight: 760;
            line-height: 1.12;
            margin: 0;
        }

        .header-subtitle {
            color: var(--muted);
            margin: 8px 0 0;
            max-width: 760px;
        }

        .header-meta {
            display: flex;
            flex-wrap: wrap;
            justify-content: flex-end;
            gap: 8px;
        }

        .chip {
            align-items: center;
            background: var(--panel);
            border: 1px solid var(--line);
            border-radius: 999px;
            color: var(--muted);
            display: inline-flex;
            font-size: 13px;
            font-weight: 650;
            gap: 7px;
            min-height: 34px;
            padding: 6px 12px;
            white-space: nowrap;
        }

        .chip strong {
            color: var(--text);
            font-weight: 760;
        }

        .toolbar {
            align-items: end;
            background: var(--panel);
            border: 1px solid var(--line);
            border-radius: 8px;
            box-shadow: var(--shadow);
            display: grid;
            gap: 16px;
            grid-template-columns: minmax(280px, 1fr);
            margin-bottom: 18px;
            padding: 16px;
        }

        .control-label {
            color: var(--muted);
            display: block;
            font-size: 12px;
            font-weight: 750;
            margin-bottom: 7px;
            text-transform: uppercase;
        }

        .btn-primary {
            align-items: center;
            background: var(--accent);
            border: 1px solid var(--accent);
            border-radius: 7px;
            color: white;
            cursor: pointer;
            display: inline-flex;
            font-weight: 720;
            gap: 8px;
            justify-content: center;
            min-height: 42px;
            padding: 9px 14px;
            transition: background-color 0.18s ease, border-color 0.18s ease, box-shadow 0.18s ease, transform 0.18s ease;
            white-space: nowrap;
        }

        .btn-primary:hover {
            background: var(--accent-dark);
            border-color: var(--accent-dark);
            box-shadow: 0 10px 24px rgba(37, 99, 235, 0.22);
            transform: translateY(-1px);
        }

        .metric-grid {
            display: grid;
            gap: 12px;
            grid-template-columns: repeat(4, minmax(0, 1fr));
            margin-bottom: 18px;
        }

        .metric-tile {
            background: var(--panel);
            border: 1px solid var(--line);
            border-radius: 8px;
            padding: 14px 16px;
        }

        .metric-label {
            color: var(--muted);
            font-size: 12px;
            font-weight: 750;
            margin-bottom: 5px;
            text-transform: uppercase;
        }

        .metric-value {
            color: var(--text);
            font-size: 22px;
            font-weight: 780;
            line-height: 1.15;
        }

        .metric-note {
            color: var(--muted);
            font-size: 12px;
            margin-top: 4px;
        }

        .panel {
            background: var(--panel);
            border: 1px solid var(--line);
            border-radius: 8px;
            box-shadow: var(--shadow);
            margin-bottom: 18px;
            overflow: hidden;
        }

        .panel-header {
            align-items: center;
            background: var(--panel);
            border-bottom: 1px solid var(--line);
            display: flex;
            flex-wrap: wrap;
            gap: 12px;
            justify-content: space-between;
            padding: 14px 16px;
        }

        .panel-header-actions {
            align-items: center;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            justify-content: flex-end;
        }

        .btn-reset {
            align-items: center;
            background: var(--panel);
            border: 1px solid var(--line-strong);
            border-radius: 7px;
            color: var(--muted);
            cursor: pointer;
            display: inline-flex;
            font-size: 12px;
            font-weight: 720;
            gap: 7px;
            justify-content: center;
            min-height: 38px;
            padding: 6px 11px;
            transition: background-color 0.18s ease, border-color 0.18s ease, color 0.18s ease, box-shadow 0.18s ease;
            white-space: nowrap;
        }

        .btn-reset:hover {
            background: var(--panel-soft);
            border-color: var(--accent);
            color: var(--accent-dark);
        }

        .btn-reset:focus-visible {
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.16);
            outline: none;
        }

        .axis-mode-toggle {
            background: var(--panel-soft);
            border: 1px solid var(--line-strong);
            border-radius: 7px;
            display: inline-flex;
            padding: 3px;
        }

        .axis-mode-button {
            background: transparent;
            border: 0;
            border-radius: 5px;
            color: var(--muted);
            cursor: pointer;
            font-size: 12px;
            font-weight: 720;
            min-height: 30px;
            padding: 5px 10px;
            transition: background-color 0.18s ease, box-shadow 0.18s ease, color 0.18s ease;
            white-space: nowrap;
        }

        .axis-mode-button:hover {
            color: var(--text);
        }

        .axis-mode-button:focus-visible {
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.2);
            outline: none;
        }

        .axis-mode-button.active {
            background: var(--panel);
            box-shadow: 0 1px 3px rgba(15, 23, 42, 0.14);
            color: var(--accent-dark);
        }

        .panel-title {
            color: var(--text);
            font-size: 15px;
            font-weight: 760;
            margin: 0;
        }

        .panel-body {
            padding: 16px;
        }

        .plot-body {
            padding: 10px 14px 16px;
        }

        #coverage-plot {
            overflow-x: auto;
            padding-bottom: 4px;
        }

        #coverage-plot svg {
            display: block;
        }

        .tooltip {
            position: absolute;
            background: rgba(15, 23, 42, 0.94);
            border: 1px solid rgba(255, 255, 255, 0.12);
            border-radius: 7px;
            box-shadow: 0 18px 36px rgba(15, 23, 42, 0.22);
            color: white;
            font-size: 12px;
            line-height: 1.45;
            max-width: 340px;
            padding: 10px 12px;
            pointer-events: none;
            z-index: 1000;
        }

        .tooltip .text-muted {
            color: #cbd5e1 !important;
        }

        .fragment,
        .fragment-indicator {
            cursor: pointer;
            transition: opacity 0.18s ease, filter 0.18s ease;
        }

        .fragment:hover,
        .fragment-indicator:hover {
            filter: brightness(1.05);
            opacity: 0.92;
        }

        .legend {
            align-items: center;
            border-top: 1px solid var(--line);
            display: flex;
            flex-wrap: wrap;
            font-size: 12px;
            gap: 8px;
            margin-top: 14px;
            padding-top: 12px;
        }

        .legend-item {
            align-items: center;
            background: var(--panel-soft);
            border: 1px solid var(--line);
            border-radius: 999px;
            color: var(--muted);
            display: inline-flex;
            gap: 7px;
            max-width: 100%;
            padding: 5px 10px;
        }

        .legend-color {
            border-radius: 999px;
            flex: 0 0 auto;
            height: 9px;
            width: 9px;
        }

        .axis line,
        .axis path {
            stroke: var(--line-strong);
        }

        .axis text {
            fill: var(--muted);
            font-size: 12px;
        }

        .brush .selection {
            fill: var(--accent);
            fill-opacity: 0.14;
            stroke: var(--accent);
            stroke-opacity: 0.72;
        }

        .annotation-item { 
            border-bottom: 1px solid var(--line);
            margin-bottom: 15px;
            padding-bottom: 10px;
        }

        .annotation-item:last-child {
            border-bottom: none;
        }

        #annotation-loading {
            padding: 24px;
            text-align: center;
        }

        .modal-content {
            border: 0;
            border-radius: 8px;
            box-shadow: 0 24px 70px rgba(15, 23, 42, 0.24);
        }

        .modal-header,
        .modal-footer {
            border-color: var(--line);
        }

        .card {
            border: 1px solid var(--line);
            border-radius: 8px;
        }

        .card-header {
            background: var(--panel-soft);
            border-bottom: 1px solid var(--line);
            color: var(--text);
            font-weight: 760;
        }

        /* Genome selector: normalize PrimeVue/Aura styling so the control and overlay
           look like one component instead of a mix of theme defaults and local CSS. */
        .genome-select.p-select {
            align-items: stretch;
            background: var(--panel) !important;
            border: 1px solid var(--line-strong) !important;
            border-radius: 8px;
            box-shadow: none;
            color: var(--text);
            min-height: 50px;
            overflow: hidden;
            transition: border-color 0.16s ease, box-shadow 0.16s ease;
            width: 100%;
        }

        .genome-select.p-select:not(.p-disabled):hover {
            border-color: #9aa9bd !important;
        }

        .genome-select.p-select.p-focus,
        .genome-select.p-select[data-p-focused="true"] {
            border-color: var(--accent) !important;
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.12) !important;
            outline: none;
        }

        .genome-select .p-select-label {
            align-items: center;
            background: transparent !important;
            color: var(--text) !important;
            display: flex;
            min-width: 0;
            padding: 8px 12px !important;
        }

        .genome-select .p-select-dropdown {
            align-items: center;
            background: var(--panel-soft) !important;
            border-left: 1px solid var(--line) !important;
            color: var(--muted) !important;
            display: flex;
            flex: 0 0 42px;
            justify-content: center;
            width: 42px;
        }

        .genome-select .p-select-dropdown:hover {
            background: #f1f5f9 !important;
            color: var(--text) !important;
        }

        .genome-select-value {
            align-items: center;
            display: flex;
            gap: 12px;
            min-width: 0;
            width: 100%;
        }

        .genome-select-value-main {
            display: flex;
            flex: 1 1 auto;
            flex-direction: column;
            gap: 2px;
            min-width: 0;
        }

        .genome-select-value-name {
            color: var(--text);
            font-size: 14px;
            font-weight: 700;
            line-height: 1.25;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }

        .genome-select-value-meta,
        .genome-option-meta {
            align-items: center;
            color: var(--muted);
            display: flex;
            flex-wrap: wrap;
            font-size: 11px;
            gap: 5px;
            line-height: 1.25;
            min-width: 0;
        }

        .genome-select-placeholder {
            color: #94a3b8;
            font-size: 14px;
            font-weight: 500;
        }

        .p-select-overlay,
        .p-select-panel {
            background: #ffffff !important;
            border: 1px solid var(--line) !important;
            border-radius: 8px !important;
            box-shadow: 0 20px 45px rgba(15, 23, 42, 0.18) !important;
            color: var(--text) !important;
            margin-top: 4px;
            overflow: hidden;
            z-index: 1100 !important;
        }

        .p-select-header,
        .p-select-filter-container {
            background: #ffffff !important;
            border-bottom: 1px solid var(--line) !important;
            padding: 10px !important;
        }

        .p-select-filter {
            background: var(--panel-soft) !important;
            border: 1px solid var(--line-strong) !important;
            border-radius: 7px !important;
            box-shadow: none !important;
            color: var(--text) !important;
            font-size: 13px;
            min-height: 38px;
            outline: none;
            padding: 8px 34px 8px 10px !important;
            width: 100%;
        }

        .p-select-filter:focus {
            background: #ffffff !important;
            border-color: var(--accent) !important;
            box-shadow: 0 0 0 3px rgba(37, 99, 235, 0.10) !important;
        }

        .p-select-list-container,
        .p-select-list {
            background: #ffffff !important;
        }

        .p-select-list {
            padding: 6px !important;
        }

        .p-select-option {
            background: #ffffff !important;
            border-radius: 6px !important;
            color: var(--text) !important;
            margin: 1px 0;
            padding: 8px 10px !important;
            transition: background-color 0.12s ease;
        }

        .p-select-option.p-focus,
        .p-select-option:hover {
            background: #f1f5f9 !important;
            color: var(--text) !important;
        }

        .p-select-option.p-select-option-selected,
        .p-select-option.p-highlight,
        .p-select-option[aria-selected="true"] {
            background: var(--accent-soft) !important;
            color: var(--text) !important;
        }

        .p-select-empty-message {
            background: #ffffff !important;
            color: var(--muted) !important;
            padding: 12px !important;
        }

        .genome-option {
            align-items: center;
            display: flex;
            gap: 12px;
            min-width: 0;
            width: 100%;
        }

        .genome-option-main {
            display: flex;
            flex: 1 1 auto;
            flex-direction: column;
            gap: 3px;
            min-width: 0;
        }

        .genome-name {
            color: var(--text);
            font-size: 13px;
            font-weight: 700;
            line-height: 1.3;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }

        .genome-taxid {
            color: var(--muted);
            white-space: nowrap;
        }

        .genome-parent {
            color: var(--teal);
            font-weight: 700;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }

        .genome-meta-separator {
            color: #cbd5e1;
        }

        .read-count {
            background: var(--accent-soft);
            border: 1px solid #bfdbfe;
            border-radius: 999px;
            color: var(--accent-dark);
            flex: 0 0 auto;
            font-size: 11px;
            font-variant-numeric: tabular-nums;
            font-weight: 760;
            padding: 3px 8px;
            white-space: nowrap;
        }

        .variant-marker {
            cursor: pointer;
            stroke-width: 0px;
            opacity: 0.5;
            transition: opacity 0.18s ease, stroke-width 0.18s ease;
        }

        .variant-marker:hover {
            opacity: 0.95;
            stroke-width: 2px;
        }

        .variant-track-label {
            fill: var(--muted);
            font-size: 11px;
            font-weight: 720;
            font-opacity: 0.8;
        }

        .variant-badge {
            background: var(--rose);
            border-radius: 999px;
            color: white;
            font-size: 12px;
            font-weight: 760;
            padding: 4px 9px;
            white-space: nowrap;
        }

        .variant-source {
            align-items: center;
            color: var(--muted);
            display: flex;
            flex-wrap: wrap;
            font-size: 13px;
            gap: 8px;
            margin-bottom: 12px;
        }

        .variant-table {
            border-color: var(--line);
            font-size: 13px;
            margin: 0;
        }

        .variant-table-wrap {
            border: 1px solid var(--line);
            border-radius: 8px;
            max-height: 430px;
            overflow: auto;
        }

        .variant-table thead th {
            background: var(--panel-soft);
            border-bottom: 1px solid var(--line);
            color: var(--muted);
            font-size: 11px;
            font-weight: 760;
            position: sticky;
            text-transform: uppercase;
            top: 0;
            white-space: nowrap;
            z-index: 1;
        }

        .variant-table td {
            color: var(--text);
            vertical-align: top;
        }

        .variant-table tbody tr {
            transition: background-color 0.16s ease;
        }

        .variant-table tbody tr:hover {
            background: #eef6ff;
        }

        .variant-row-selected {
            outline: 2px solid var(--accent);
            outline-offset: -2px;
        }

        .detail-grid {
            display: grid;
            gap: 12px;
            grid-template-columns: repeat(2, minmax(0, 1fr));
        }

        .detail-item {
            background: var(--panel-soft);
            border: 1px solid var(--line);
            border-radius: 7px;
            padding: 12px;
        }

        .detail-label {
            color: var(--muted);
            font-size: 11px;
            font-weight: 760;
            text-transform: uppercase;
        }

        .detail-value {
            color: var(--text);
            font-weight: 720;
            margin-top: 4px;
            overflow-wrap: anywhere;
        }

        .profile-subtitle {
            color: var(--muted);
            font-size: 13px;
            line-height: 1.5;
            margin: 6px 0 0;
        }

        .profile-panel-header {
            align-items: flex-start;
            flex-wrap: wrap;
        }

        .profile-panel-header > div {
            flex: 1 1 280px;
        }

        .profile-status {
            background: #eef2f6;
            border: 1px solid var(--line);
            border-radius: 999px;
            color: #475569;
            font-size: 12px;
            font-weight: 700;
            padding: 5px 12px;
            white-space: nowrap;
        }

        .profile-status.qualified {
            background: #ecfdf5;
            border-color: #a7f3d0;
            color: #065f46;
        }

        .profile-status.filtered {
            background: #fff7ed;
            border-color: #fed7aa;
            color: #9a3412;
        }

        .profile-highlights {
            display: grid;
            gap: 12px;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            margin-bottom: 18px;
        }

        .profile-highlight .detail-value {
            font-size: 25px;
            font-variant-numeric: tabular-nums;
        }

        .profile-content-grid {
            display: grid;
            gap: 28px;
            grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
            align-items: start;
        }

        .profile-groups {
            display: flex;
            flex-direction: column;
            gap: 22px;
            min-width: 0;
            order: 2;
        }

        .profile-network-column {
            border-right: 1px solid var(--line);
            min-width: 0;
            order: 1;
            padding-right: 28px;
        }

        .profile-network-header {
            align-items: flex-start;
            display: flex;
            gap: 12px;
            justify-content: space-between;
            margin-bottom: 10px;
        }

        .profile-group h3,
        .profile-notes h3 {
            color: var(--text);
            font-size: 14px;
            font-weight: 750;
            margin: 0 0 10px;
        }

        .profile-group dl {
            margin: 0;
        }

        .profile-row {
            align-items: baseline;
            border-bottom: 1px solid var(--line);
            display: grid;
            gap: 12px;
            grid-template-columns: minmax(0, 1fr) minmax(0, 1fr);
            padding: 9px 0;
        }

        .profile-row dt {
            color: var(--muted);
            font-size: 13px;
            font-weight: 500;
        }

        .profile-row dd {
            color: var(--text);
            font-size: 13px;
            font-variant-numeric: tabular-nums;
            font-weight: 650;
            margin: 0;
            overflow-wrap: anywhere;
            text-align: right;
        }

        .profile-notes {
            background: var(--panel-soft);
            border-left: 3px solid var(--line-strong);
            border-radius: 4px;
            margin-top: 20px;
            padding: 14px 16px;
        }

        .profile-notes p {
            font-size: 13px;
            line-height: 1.6;
            margin: 0;
            overflow-wrap: anywhere;
            white-space: pre-wrap;
        }

        .empty-state {
            align-items: center;
            background: var(--panel-soft);
            border: 1px dashed var(--line-strong);
            border-radius: 8px;
            color: var(--muted);
            display: flex;
            min-height: 92px;
            padding: 18px;
        }

        .empty-state code {
            background: var(--accent-soft);
            border-radius: 5px;
            color: var(--accent-dark);
            padding: 2px 5px;
        }

        .network-subtitle {
            color: var(--muted);
            font-size: 13px;
            line-height: 1.5;
            margin: 5px 0 0;
        }

        #network-plot {
            min-height: 520px;
            position: relative;
            width: 100%;
        }

        #network-plot svg {
            display: block;
            min-height: 520px;
            width: 100%;
        }

        .network-link {
            fill: none;
            stroke: #8a98aa;
            stroke-opacity: 0.46;
        }

        .network-node {
            cursor: grab;
            stroke: #ffffff;
            stroke-width: 2px;
        }

        .network-node.selected {
            stroke: #172033;
            stroke-width: 3px;
        }

        .network-node.clickable {
            cursor: pointer;
        }

        .network-label {
            fill: var(--text);
            font-size: 11px;
            paint-order: stroke;
            pointer-events: none;
            stroke: rgba(255, 255, 255, 0.94);
            stroke-width: 4px;
            stroke-linejoin: round;
        }

        .network-label.selected {
            font-weight: 800;
        }

        .network-legend {
            align-items: center;
            border-top: 1px solid var(--line);
            display: flex;
            flex-wrap: wrap;
            font-size: 12px;
            gap: 8px;
            margin-top: 10px;
            padding-top: 12px;
        }

        .network-species-panel {
            margin-top: 18px;
        }

        .network-species-title {
            color: var(--text);
            font-size: 13px;
            font-weight: 750;
            margin: 0 0 8px;
        }

        .network-species-table {
            border-color: var(--line);
            font-size: 13px;
            margin: 0;
        }

        .network-species-table-wrap {
            border: 1px solid var(--line);
            border-radius: 8px;
            max-height: 260px;
            overflow: auto;
        }

        .network-species-table thead th {
            background: var(--panel-soft);
            border-bottom: 1px solid var(--line);
            color: var(--muted);
            font-size: 11px;
            font-weight: 760;
            position: sticky;
            text-transform: uppercase;
            top: 0;
            white-space: nowrap;
            z-index: 1;
        }

        .network-species-table td {
            color: var(--text);
            vertical-align: top;
        }

        .network-species-table tbody tr {
            transition: background-color 0.16s ease;
        }

        .network-species-table tbody tr:hover {
            background: #eef6ff;
        }

        .network-species-table tr.network-species-row-selected {
            background: var(--accent-soft);
        }

        @media (max-width: 900px) {
            .app-shell {
                padding: 18px;
            }

            .app-header {
                align-items: flex-start;
                flex-direction: column;
            }

            .header-meta {
                justify-content: flex-start;
            }

            .toolbar,
            .metric-grid,
            .detail-grid,
            .profile-content-grid {
                grid-template-columns: 1fr;
            }

            .profile-network-column {
                border-right: 0;
                border-top: 1px solid var(--line);
                padding-right: 0;
                padding-top: 24px;
            }

            .profile-highlights {
                grid-template-columns: repeat(2, minmax(0, 1fr));
            }

            .btn-primary {
                width: 100%;
            }

            .genome-select-value,
            .genome-option {
                gap: 8px;
            }

            .genome-select-value .read-count,
            .genome-option .read-count {
                font-size: 10px;
                padding: 2px 6px;
            }
        }
    </style>
</head>
<body>
    <div id="page-loader">
        <div class="page-loader-spinner"></div>
        <div class="page-loader-text">Loading report...</div>
    </div>
    <div id="app" class="app-shell">
        <header class="app-header">
            <div>
                <p class="eyebrow">SPADES-GOTTCHA2</p>
                <h1>Genome Signature Browser</h1>
                <p class="header-subtitle" v-if="currentGenome">
                    {{ currentGenome.name }}<span v-if="currentGenome.parentName"> ({{ currentGenome.parentName }})</span> - {{ currentGenome.taxid }}
                </p>
            </div>
            <div class="header-meta">
                <span class="chip"><i class="pi pi-database"></i><strong>{{ genomeOptions.length.toLocaleString() }}</strong> genomes</span>
                <span class="chip"><i class="pi pi-file"></i><strong>{{ vcfFiles.length.toLocaleString() }}</strong> VCF sources</span>
            </div>
        </header>

        <section class="toolbar">
            <div>
                <label class="control-label">Genome</label>
                <p-select
                    v-model="selectedGenome"
                    :options="genomeOptions"
                    option-label="label"
                    option-value="value"
                    placeholder="Search or select a genome..."
                    filter
                    filter-placeholder="Search name, parent, or TaxID"
                    scroll-height="24rem"
                    @change="handleGenomeSelect"
                    class="genome-select"
                >
                    <template #value>
                        <div v-if="selectedGenomeOption" class="genome-select-value">
                            <div class="genome-select-value-main">
                                <span class="genome-select-value-name">{{ selectedGenomeOption.name }}</span>
                                <span class="genome-select-value-meta">
                                    <span v-if="selectedGenomeOption.parentName" class="genome-parent">{{ selectedGenomeOption.parentName }}</span>
                                    <span v-if="selectedGenomeOption.parentName" class="genome-meta-separator">&bull;</span>
                                    <span class="genome-taxid">TaxID {{ selectedGenomeOption.taxid }}</span>
                                </span>
                            </div>
                            <span class="read-count">{{ selectedGenomeOption.readCount.toLocaleString() }} reads</span>
                        </div>
                        <span v-else class="genome-select-placeholder">Search or select a genome...</span>
                    </template>
                    <template #option="slotProps">
                        <div class="genome-option">
                            <div class="genome-option-main">
                                <span class="genome-name">{{ slotProps.option.name }}</span>
                                <span class="genome-option-meta">
                                    <span v-if="slotProps.option.parentName" class="genome-parent">{{ slotProps.option.parentName }}</span>
                                    <span v-if="slotProps.option.parentName" class="genome-meta-separator">&bull;</span>
                                    <span class="genome-taxid">TaxID {{ slotProps.option.taxid }}</span>
                                </span>
                            </div>
                            <span class="read-count">{{ slotProps.option.readCount.toLocaleString() }} reads</span>
                        </div>
                    </template>
                </p-select>
            </div>
        </section>

        <section class="metric-grid" v-if="currentSummary">
            <div class="metric-tile">
                <div class="metric-label">Reads</div>
                <div class="metric-value">{{ currentSummary.readCount.toLocaleString() }}</div>
                <div class="metric-note">mapped to signatures</div>
            </div>
            <div class="metric-tile">
                <div class="metric-label">Fragments</div>
                <div class="metric-value">{{ currentSummary.fragmentCount.toLocaleString() }}</div>
                <div class="metric-note">shown in this genome</div>
            </div>
            <div class="metric-tile">
                <div class="metric-label">Signature Coverage</div>
                <div class="metric-value">{{ currentSummary.overallCoverage.toFixed(1) }}%</div>
                <div class="metric-note">{{ currentSummary.coverageSource }}</div>
            </div>
            <div class="metric-tile">
                <div class="metric-label">Mean Depth</div>
                <div class="metric-value">{{ currentSummary.meanDepth.toFixed(1) }}x</div>
                <div class="metric-note">across fragments</div>
            </div>
        </section>
        
        <section class="panel">
            <div class="panel-header">
                <h2 class="panel-title">Signature Coverage</h2>
                <div class="panel-header-actions">
                    <span class="chip" v-if="currentGenome"><i class="pi pi-chart-bar"></i> <strong>{{ coverageAxisLength.toLocaleString() }} bp</strong></span>
                    <div class="axis-mode-toggle" role="group" aria-label="Signature coverage x-axis scale">
                        <button
                            type="button"
                            class="axis-mode-button"
                            :class="{ active: coverageAxisMode === 'signature' }"
                            :aria-pressed="coverageAxisMode === 'signature'"
                            @click="selectCoverageAxisMode('signature')"
                        >
                            Total Signatures
                        </button>
                        <button
                            type="button"
                            class="axis-mode-button"
                            :class="{ active: coverageAxisMode === 'genome' }"
                            :aria-pressed="coverageAxisMode === 'genome'"
                            @click="selectCoverageAxisMode('genome')"
                        >
                            Genomic Position
                        </button>
                    </div>
                    <button @click="resetZoom" class="btn-reset" type="button" title="Reset coverage zoom">
                        <i class="pi pi-refresh"></i>
                        Reset
                    </button>
                </div>
            </div>
            <div class="panel-body plot-body">
                <div id="coverage-plot"></div>
                <div class="legend" id="legend"></div>
            </div>
        </section>

        <section class="panel" id="variant-card">
            <div class="panel-header">
                <h2 class="panel-title">VCF Variants</h2>
                <span id="variant-count-badge" class="variant-badge">0 variants</span>
            </div>
            <div class="panel-body">
                <div id="variant-file-list" class="variant-source"></div>
                <div id="variant-table"></div>
            </div>
        </section>
        
        <section class="panel" v-if="currentGenome">
            <div class="panel-header profile-panel-header">
                <div>
                    <h2 class="panel-title">Profiling Results</h2>
                    <p class="profile-subtitle">Full profiling results for this genome. Read counts and depth may differ from the alignment summary above. An em dash (—) means the value is unavailable.</p>
                </div>
                <span class="profile-status" :class="profileStatus.tone">{{ profileStatus.label }}</span>
            </div>
            <div class="panel-body">
                <div class="profile-content-grid">
                    <div>
                        <div class="profile-groups">
                            <section class="profile-group" v-for="group in profileSections" :key="group.title">
                                <h3>{{ group.title }}</h3>
                                <dl>
                                    <div class="profile-row" v-for="metric in group.metrics" :key="metric.field" :title="metric.hint">
                                        <dt>{{ metric.label }}</dt>
                                        <dd>{{ formatProfileValue(currentGenome.profile[metric.field], metric.format) }}</dd>
                                    </div>
                                </dl>
                            </section>
                        </div>
                        <div class="profile-notes" v-if="currentGenome.profile?.NOTE">
                            <h3>Profiling notes</h3>
                            <p>{{ currentGenome.profile.NOTE }}</p>
                        </div>
                    </div>
                    <aside class="profile-network-column">
                        <div class="profile-highlights">
                            <div class="detail-item profile-highlight" v-for="metric in profileHighlights" :key="metric.field" :title="metric.hint">
                                <div class="detail-label">{{ metric.label }}</div>
                                <div class="detail-value">{{ formatProfileValue(currentGenome.profile?.[metric.field], metric.format) }}</div>
                            </div>
                        </div>
                        <div class="profile-network-header">
                            <div>
                                <h3 class="panel-title">Cross-species Network</h3>
                                <p class="network-subtitle">Primary → other cross-TaxID alignments. Nodes from the same species share a color.</p>
                            </div>
                            <span id="network-count-badge" class="profile-status">0 nodes · 0 edges</span>
                        </div>
                        <div id="network-plot"></div>
                        <div class="network-legend" id="network-legend"></div>
                        <div class="network-species-panel">
                            <h4 class="network-species-title">Relevant species profiled</h4>
                            <div id="network-species-table"></div>
                        </div>
                    </aside>
                </div>
            </div>
        </section>
        <section class="panel" v-else>
            <div class="panel-header">
                <h2 class="panel-title">Profiling Results</h2>
            </div>
            <div class="panel-body">
                <div class="empty-state">Select a genome to view information</div>
            </div>
        </section>
    </div>

    <!-- NCBI Annotation Query Modal -->
    <div class="modal fade" id="annotationModal" tabindex="-1" aria-labelledby="annotationModalLabel" aria-hidden="true">
        <div class="modal-dialog modal-lg">
            <div class="modal-content">
                <div class="modal-header">
                    <h5 class="modal-title" id="annotationModalLabel">NCBI Nucleotide Annotation Query</h5>
                    <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                </div>
                <div class="modal-body">
                    <div id="annotation-loading" class="d-none">
                        <div class="spinner-border" role="status">
                            <span class="visually-hidden">Loading...</span>
                        </div>
                        <p>Fetching annotations...</p>
                    </div>
                    <div id="annotation-results"></div>
                </div>
                <div class="modal-footer">
                    <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
                </div>
            </div>
        </div>
    </div>

    <!-- Bootstrap 5 JS Bundle with Popper -->
    <script src="/publicdata/js/bootstrap.bundle.min.js"></script>
    <script>
    window.browserResourcesReady.then(function() {
        // The data will be injected here
        const genomeData = GENOME_DATA_PLACEHOLDER;
        const speciesData = SPECIES_DATA_PLACEHOLDER;
        const coverageData = COVERAGE_DATA_PLACEHOLDER;
        const variantData = VARIANT_DATA_PLACEHOLDER;
        const vcfFiles = VCF_FILE_DATA_PLACEHOLDER;
        const networkNodes = NETWORK_NODE_DATA_PLACEHOLDER;
        const networkEdges = NETWORK_EDGE_DATA_PLACEHOLDER;
        const minVariantDepth = MIN_VARIANT_DEPTH_PLACEHOLDER;
        const variantsBySeq = variantData.reduce((acc, variant) => {
            const seqName = String(variant.seq_name || '');
            if (!seqName) return acc;
            if (!acc[seqName]) acc[seqName] = [];
            acc[seqName].push(variant);
            return acc;
        }, {});
        
        // Initialize D3 variables
        let svg, xScale, yScale, xAxis, yAxis, xAxisTop, brush;
        let coveragePlotObserver, coverageResizeFrame;
        const margin = {top: 50, right: 30, bottom: 50, left: 60};
        const minimumChartWidth = 640;
        let width = 1100 - margin.left - margin.right;
        const height = 400 - margin.top - margin.bottom;
        
        // Initialize Vue app
        const { createApp, ref, computed, onMounted, watch } = Vue;
        
        // Bootstrap modal instance
        let annotationModal;
        
        function initializeAnnotationModal() {
            annotationModal = new bootstrap.Modal(document.getElementById('annotationModal'));
        }
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', initializeAnnotationModal);
        } else {
            initializeAnnotationModal();
        }
        
        // Reactive state variables accessible globally
        let selectedTaxid, selectedGenome, currentGenome, currentFragment, coverageAxisMode;
        let networkSimulation = null;

        const profileHighlights = [
            { field: 'READ_COUNT', label: 'Profiled reads', format: 'integer', hint: 'Read count from the full profiling report.' },
            { field: 'SNI_SCORE', label: 'SNI score', format: 'score', hint: 'Signature nucleotide identity score used for taxonomic filtering.' },
        ];
        const profileGroups = [
            { title: 'Taxonomy & reference', metrics: [
                { field: 'NAME', label: 'Name' },
                { field: 'TAXID', label: 'Taxonomy ID' },
                { field: 'PARENT_NAME', label: 'Parent taxon' },
                { field: 'PARENT_TAXID', label: 'Parent taxonomy ID' },
                { field: 'SIG_LEVEL', label: 'Signature level' },
                { field: 'GENOME_SIZE', label: 'Genome size', format: 'bases' },
                { field: 'GENOME_COUNT', label: 'Reference entries', format: 'integer' },
            ]},
            { title: 'Coverage & depth', metrics: [
                { field: 'TOTAL_SIG_LEN', label: 'Total signature length', format: 'bases' },
                { field: 'COVERED_SIG_LEN', label: 'Covered signature length', format: 'bases' },
                { field: 'MAPPED_SIG_LEN', label: 'Mapped signature length', format: 'bases' },
                { field: 'SIG_COV_RAW', label: 'Aggregate signature coverage', format: 'percent', hint: 'Covered signature length divided by total signature length.' },
                { field: 'DEPTH', label: 'Profile depth', format: 'depth' },
                { field: 'COVERED_SIG_DEPTH', label: 'Depth on covered signatures', format: 'depth' },
                { field: 'COVERED_MAPPED_SIG_COV', label: 'Coverage of mapped signatures', format: 'percent' },
                { field: 'ZSCORE', label: 'Depth z-score', format: 'number' },
            ]},
            { title: 'Identity & read evidence', metrics: [
                { field: 'SNI_CI95_LH', label: 'SNI 95% confidence interval' },
                { field: 'ALN_IDENTITY', label: 'Alignment identity', format: 'percent' },
                { field: 'CONSENSUS_SEQ_IDENTITY', label: 'Consensus identity', format: 'percent', hint: 'Consensus identity before sequencing-error correction.' },
                { field: 'TOTAL_BP_MAPPED', label: 'Mapped bases', format: 'bases' },
                { field: 'TOTAL_READ_LEN', label: 'Aligned query length', format: 'bases' },
                { field: 'TOTAL_BP_MISMATCH', label: 'Mismatched bases', format: 'bases' },
                { field: 'TOTAL_BP_INDEL', label: 'Inserted / deleted bases', format: 'bases' },
                { field: 'AOI_READ_COUNT', label: 'Reads matching signatures of interest', format: 'integer' },
            ]},
            { title: 'Abundance estimates', metrics: [
                { field: 'REL_ABUNDANCE_DEPTH', label: 'Relative abundance by depth', format: 'percent' },
                { field: 'REL_ABUNDANCE_GC', label: 'Relative abundance by genomic content', format: 'percent' },
            ]},
        ];

        function formatProfileValue(value, format) {
            if (value === null || value === undefined || value === '') return '—';
            if (!format) return String(value);
            const numeric = Number(value);
            if (!Number.isFinite(numeric)) return '—';
            const number = format === 'percent' ? numeric * 100 : numeric;
            const suffix = { percent: '%', bases: ' bp', depth: '×' }[format] || '';
            // Preserve small nonzero abundances instead of rounding them to zero.
            if (number !== 0 && Math.abs(number) < 0.0001) {
                return `${number.toExponential(2)}${suffix}`;
            }
            const digits = ['integer', 'bases'].includes(format) ? 0 : format === 'score' ? 6 : 4;
            return number.toLocaleString(undefined, { maximumFractionDigits: digits }) + suffix;
        }

        function getProfileStatus(profile) {
            if (!Object.prototype.hasOwnProperty.call(profile, 'NOTE')) {
                return { label: 'Status unavailable', tone: 'unknown' };
            }
            const note = profile.NOTE || '';
            if (note.includes('species SNI_SCORE')) return { label: 'Filtered out', tone: 'filtered' };
            return { label: 'Qualified', tone: 'qualified' };
        }
        
        // Register Vue-Select component
        const app = createApp({
            setup() {
                // Reactive state - assign to outer variables for global access
                selectedTaxid = ref('');
                currentGenome = ref(null);
                currentFragment = ref(null);
                coverageAxisMode = ref('signature');
                selectedGenome = ref(null);
                
                // Calculate unique taxids and read counts
                const uniqueTaxidsInCoverage = [...new Set(coverageData.map(item => item.genome_taxid))];
                
                // Calculate total read counts for each genome
                const genomeTotalReads = {};
                coverageData.forEach(item => {
                    if (!genomeTotalReads[item.genome_taxid]) {
                        genomeTotalReads[item.genome_taxid] = 0;
                    }
                    genomeTotalReads[item.genome_taxid] += item.numreads;
                });
                
                // Prepare genomes with read counts
                const genomesWithReadCounts = genomeData
                    .filter(genome => uniqueTaxidsInCoverage.includes(genome.taxid))
                    .map(genome => ({
                        ...genome,
                        readCount: genomeTotalReads[genome.taxid] || 0
                    }))
                    .sort((a, b) => b.readCount - a.readCount); // Sort by read count descending
                
                // Format options for vue-select
                const genomeOptions = genomesWithReadCounts.map(genome => ({
                    value: genome.taxid,
                    label: genome.parentName
                        ? `${genome.name} (${genome.parentName}) (${genome.taxid})`
                        : `${genome.name} (${genome.taxid})`,
                    name: genome.name,
                    parentName: genome.parentName || '',
                    taxid: genome.taxid,
                    readCount: genome.readCount
                }));

                const selectedGenomeOption = computed(() => {
                    const value = selectedGenome.value;
                    if (value === null || value === undefined || value === '') return null;
                    return genomeOptions.find(option => String(option.value) === String(value)) || null;
                });

                function percentValue(value) {
                    if (value === null || value === undefined || value === '') return null;
                    const numeric = Number(value);
                    if (!Number.isFinite(numeric)) return null;
                    return numeric <= 1 ? numeric * 100 : numeric;
                }

                function formatPercent(value) {
                    const pct = percentValue(value);
                    return pct === null ? 'N/A' : `${pct.toFixed(1)}%`;
                }

                const currentSummary = computed(() => {
                    const genome = currentGenome.value;
                    if (!genome) return null;

                    const fragments = coverageData.filter(d => String(d.genome_taxid) === String(genome.taxid));
                    const fragmentCount = fragments.length;
                    const readCount = genomeTotalReads[genome.taxid] || 0;
                    const coveredBases = fragments.reduce((sum, fragment) => {
                        const length = Math.max(0, Number(fragment.end_position) - Number(fragment.start_position) + 1);
                        return sum + (length * Number(fragment.coverage || 0) / 100);
                    }, 0);
                    const totalBases = fragments.reduce((sum, fragment) => {
                        return sum + Math.max(0, Number(fragment.end_position) - Number(fragment.start_position) + 1);
                    }, 0);
                    const profileCoverage = percentValue(genome.sigCov);
                    const overallCoverage = profileCoverage !== null
                        ? profileCoverage
                        : (totalBases ? (coveredBases / totalBases) * 100 : 0);
                    const coverageSource = profileCoverage !== null ? 'overall coverage' : 'weighted signatures';
                    const meanDepth = fragmentCount
                        ? fragments.reduce((sum, fragment) => sum + Number(fragment.meandepth || 0), 0) / fragmentCount
                        : 0;

                    return { fragmentCount, readCount, overallCoverage, coverageSource, meanDepth };
                });

                const coverageAxisLength = computed(() => {
                    if (!currentGenome.value) return 0;
                    return getFullXAxisDomain(currentGenome.value)[1];
                });

                const profileStatus = computed(() => getProfileStatus(currentGenome.value?.profile || {}));
                const profileSections = computed(() => {
                    const profile = currentGenome.value?.profile || {};
                    return profileGroups.map(group => ({
                        ...group,
                        metrics: group.metrics.filter(metric => Object.prototype.hasOwnProperty.call(profile, metric.field))
                    })).filter(group => group.metrics.length);
                });
                
                // Methods
                function handleGenomeChange() {
                    console.log("Selected TaxID:", selectedTaxid.value);
                    if (selectedTaxid.value) {
                        currentGenome.value = genomeData.find(g => g.taxid === selectedTaxid.value);
                        updateVisualization();
                    }
                }

                function handleGenomeSelect(eventOrValue) {
                    const value = eventOrValue && typeof eventOrValue === "object" && "value" in eventOrValue
                        ? eventOrValue.value
                        : eventOrValue;

                    if (!value) return;

                    selectedTaxid.value = String(value);
                    selectedGenome.value = selectedTaxid.value;
                    currentGenome.value = genomeData.find(g => String(g.taxid) === selectedTaxid.value);
                    updateVisualization();
                }

                function resetZoom() {
                    if (currentGenome.value) {
                        xScale.domain(getFullXAxisDomain(currentGenome.value));
                        svg.select(".brush").call(brush.move, null);
                        updateChart();
                    }
                }

                function selectCoverageAxisMode(mode) {
                    changeCoverageAxisMode(mode);
                }
                
                onMounted(() => {
                    initVisualization();
                    if (genomeOptions.length > 0) {
                        selectedGenome.value = genomeOptions[0].value;
                        selectedTaxid.value = genomeOptions[0].value;
                        currentGenome.value = genomeData.find(g => g.taxid === selectedTaxid.value);
                        updateVisualization();
                    }
                    initializeCoveragePlotObserver();
                });

                return {
                    selectedTaxid,
                    selectedGenome,
                    selectedGenomeOption,
                    genomeOptions,
                    currentGenome,
                    currentFragment,
                    currentSummary,
                    coverageAxisMode,
                    coverageAxisLength,
                    genomesWithReadCounts,
                    vcfFiles,
                    minVariantDepth,
                    formatPercent,
                    profileHighlights,
                    profileSections,
                    profileStatus,
                    formatProfileValue,
                    handleGenomeChange,
                    handleGenomeSelect,
                    selectCoverageAxisMode,
                    resetZoom
                };
            },
            components: {
                'p-select': PrimeVue.Select,
            }
        });
        
        // Configure PrimeVue with theme
        app.use(PrimeVue.Config, {
            theme: {
                preset: PrimeUIX.Themes.Aura,
                options: {
                    darkModeSelector: '.dark-mode'
                }
            }
        });
        
        // Mount Vue app
        app.mount('#app');

        // Reveal the rendered report now that Vue has mounted, and drop the loading overlay.
        document.getElementById('app').style.visibility = 'visible';
        const pageLoader = document.getElementById('page-loader');
        if (pageLoader) pageLoader.remove();
        
        // Create tooltip
        const tooltip = d3.select('body').append('div')
            .attr('class', 'tooltip')
            .style('opacity', 0);

        const chartTransitionDuration = 750;

        function usesGenomeXAxis() {
            return coverageAxisMode && coverageAxisMode.value === 'genome';
        }

        function getFullXAxisDomain(genome) {
            if (!genome) return [0, 1];

            let axisLength = Number(usesGenomeXAxis() ? genome.genomeSize : genome.totalLength);
            if (!Number.isFinite(axisLength) || axisLength <= 0) {
                const fragments = coverageData.filter(d => String(d.genome_taxid) === String(genome.taxid));
                axisLength = usesGenomeXAxis()
                    ? fragments.reduce((maximum, fragment) => Math.max(maximum, Number(fragment.end_position) || 0), 0)
                    : fragments.reduce((total, fragment) => {
                        return total + Math.max(0, Number(fragment.end_position) - Number(fragment.start_position) + 1);
                    }, 0);
            }

            return [0, Math.max(1, axisLength)];
        }

        function getFragmentXStart(fragment) {
            // Coverage coordinates are 1-based and inclusive; use half-open plot bounds.
            return usesGenomeXAxis() ? Math.max(0, Number(fragment.start_position) - 1) : Number(fragment.x_start);
        }

        function getFragmentXEnd(fragment) {
            return usesGenomeXAxis() ? Number(fragment.end_position) : Number(fragment.x_end);
        }

        function getFragmentPixelWidth(fragment) {
            return Math.max(2, xScale(getFragmentXEnd(fragment)) - xScale(getFragmentXStart(fragment)));
        }

        function getVariantXPosition(variant) {
            return usesGenomeXAxis() ? Math.max(0, Number(variant.genome_x) - 1) : Number(variant.signature_x);
        }

        function getBottomXAxisLabel() {
            return usesGenomeXAxis() ? 'Genomic Position (bp)' : 'Total Signature (bp)';
        }

        function getTopXAxisLabel() {
            return usesGenomeXAxis() ? 'Genome Position (%)' : 'Covered Signature Fragment (%)';
        }

        function changeCoverageAxisMode(mode) {
            if (!coverageAxisMode || (mode !== 'signature' && mode !== 'genome')) return;
            if (coverageAxisMode.value === mode) return;

            coverageAxisMode.value = mode;
            if (!currentGenome || !currentGenome.value || !svg || !xScale) return;

            xScale.domain(getFullXAxisDomain(currentGenome.value));
            svg.select('.brush').call(brush.move, null);
            updateChart();
        }

        function getCoveragePlotOuterWidth() {
            const plot = document.getElementById('coverage-plot');
            const availableWidth = plot ? Math.floor(plot.getBoundingClientRect().width) : 0;
            return Math.max(minimumChartWidth, availableWidth);
        }

        function initializeCoveragePlotObserver() {
            const plot = document.getElementById('coverage-plot');
            if (!plot || coveragePlotObserver || typeof ResizeObserver === 'undefined') return;

            let observedWidth = Math.floor(plot.getBoundingClientRect().width);
            coveragePlotObserver = new ResizeObserver(entries => {
                const nextWidth = Math.floor(entries[0].contentRect.width);
                if (!nextWidth || Math.abs(nextWidth - observedWidth) < 2) return;
                observedWidth = nextWidth;

                if (coverageResizeFrame) cancelAnimationFrame(coverageResizeFrame);
                coverageResizeFrame = requestAnimationFrame(() => {
                    if (currentGenome && currentGenome.value) {
                        updateVisualization();
                    } else {
                        initVisualization();
                    }
                });
            });
            coveragePlotObserver.observe(plot);
        }
        
        // Initialize the visualization container
        function initVisualization() {
            d3.select("#coverage-plot").html("");

            const chartOuterWidth = getCoveragePlotOuterWidth();
            width = chartOuterWidth - margin.left - margin.right;
            
            svg = d3.select("#coverage-plot")
                .append("svg")
                .attr("viewBox", `0 0 ${chartOuterWidth} ${height + margin.top + margin.bottom}`)
                .attr("width", chartOuterWidth)
                .attr("height", height + margin.top + margin.bottom)
                .attr("role", "img")
                .append("g")
                .attr("transform", `translate(${margin.left},${margin.top})`);
            
            // Create scales
            xScale = d3.scaleLinear().range([0, width]);
            yScale = d3.scaleLinear()
                .domain([0, 0.1]) // Will be adjusted based on data
                .range([height, 0]);
            
            // Create axes
            xAxis = svg.append("g")
                .attr("class", "axis x-axis")
                .attr("transform", `translate(0,${height})`);
            
            xAxisTop = svg.append("g")
                .attr("class", "axis x-axis-top")
                .attr("transform", "translate(0,0)");
            
            yAxis = svg.append("g")
                .attr("class", "axis y-axis");
            
            // Add axis labels
            svg.append("text")
                .attr("class", "x-axis-label")
                .attr("transform", `translate(${width/2},${height + 40})`)
                .style("text-anchor", "middle")
                .style("font-size", "12px")
                .style("font-weight", "700")
                .style("fill", "#697386")
                .text(getBottomXAxisLabel());
            
            svg.append("text")
                .attr("class", "x-axis-top-label")
                .attr("transform", `translate(${width/2},-30)`)
                .style("text-anchor", "middle")
                .style("font-size", "12px")
                .style("font-weight", "700")
                .style("fill", "#697386")
                .text(getTopXAxisLabel());
            
            svg.append("text")
                .attr("transform", "rotate(-90)")
                .attr("y", 0 - margin.left)
                .attr("x", 0 - (height / 2) - 5)
                .attr("dy", "1em")
                .style("fill", "#697386")
                .style("text-anchor", "middle")
                .style("font-size", "12px")
                .style("font-weight", "700")
                .text("Mean Depth (x)");
            
            // Create brush for zooming
            brush = d3.brushX()
                .extent([[0, 0], [width, height]])
                .on("end", brushended);
            
            svg.append("g")
                .attr("class", "brush")
                .call(brush);
        }

        function escapeHtml(value) {
            if (value === null || value === undefined) return '';
            return String(value)
                .replace(/&/g, '&amp;')
                .replace(/</g, '&lt;')
                .replace(/>/g, '&gt;')
                .replace(/"/g, '&quot;')
                .replace(/'/g, '&#039;');
        }

        function formatNullableNumber(value) {
            return value === null || value === undefined || value === '' ? '.' : Number(value).toLocaleString();
        }

        function formatNullablePercent(value) {
            if (value === null || value === undefined || value === '') return '.';
            const pct = Number(value);
            return Number.isFinite(pct) ? `${pct.toFixed(1)}%` : '.';
        }

        function formatVariantTooltip(d) {
            return `
                <strong>Variant:</strong> ${escapeHtml(d.type || 'VAR')}<br>
                <strong>Sequence:</strong> ${escapeHtml(d.seq_name)}<br>
                <strong>Position:</strong> ${Number(d.pos).toLocaleString()}${d.end && d.end !== d.pos ? '-' + Number(d.end).toLocaleString() : ''}<br>
                <strong>REF:</strong> ${escapeHtml(d.ref)}<br>
                <strong>ALT:</strong> ${escapeHtml(d.alt)}<br>
                <strong>DEPTH:</strong> ${formatNullableNumber(d.depth)}<br>
                <strong>REF_DP:</strong> ${formatNullablePercent(d.ref_dp)}<br>
                <strong>ALT_DP:</strong> ${formatNullablePercent(d.alt_dp)}<br>
                <small class="text-muted">Click to highlight this variant in the table</small>
            `;
        }

        function getMappedVariants(genomeFragments) {
            const mapped = [];
            const seen = new Set();

            genomeFragments.forEach(fragment => {
                const seqName = String(fragment.seq_name);
                const fragmentVariants = variantsBySeq[seqName] || [];
                const fragmentStart = Number(fragment.start_position);
                const fragmentEnd = Number(fragment.end_position);

                fragmentVariants.forEach(variant => {
                    const variantStart = Number(variant.pos);
                    const variantEnd = Number(variant.end || variant.pos);

                    if (variantEnd < fragmentStart || variantStart > fragmentEnd) return;

                    const markerPosition = Math.min(Math.max(variantStart, fragmentStart), fragmentEnd);
                    const signatureX = Number(fragment.x_start) + (markerPosition - fragmentStart);
                    const key = `${variant.source_index || 0}|${seqName}|${variant.pos}|${variant.end}|${variant.ref}|${variant.alt}|${fragment.x_start}`;
                    if (seen.has(key)) return;
                    seen.add(key);

                    mapped.push({
                        ...variant,
                        fragment_start: fragmentStart,
                        fragment_end: fragmentEnd,
                        relative_position: markerPosition - fragmentStart + 1,
                        signature_x: signatureX,
                        genome_x: markerPosition
                    });
                });
            });

            mapped.sort((a, b) => a.signature_x - b.signature_x || Number(a.pos) - Number(b.pos));
            mapped.forEach((variant, index) => {
                variant.variant_index = index;
            });
            return mapped;
        }

        function drawVariants(mappedVariants) {
            const variantTrackY = 24;

            if (!mappedVariants.length) return;

            svg.append('text')
                .attr('class', 'variant-track-label')
                .attr('x', 0)
                .attr('y', variantTrackY + 5)
                .text('Variants');

            svg.append('text')
                .attr('class', 'variant-track-label')
                .attr('x', 0)
                .attr('y', variantTrackY - 12)
                .text('Signature Fragments');

            const variantTypes = [...new Set(mappedVariants.map(d => d.type || 'VAR'))];
            const variantColorScale = d3.scaleOrdinal(d3.schemeSet2).domain(variantTypes);
            const markerSymbol = d3.symbol().type(d3.symbolTriangle).size(70);

            svg.selectAll('.variant-marker')
                .data(mappedVariants)
                .enter()
                .append('path')
                .attr('class', 'variant-marker')
                .attr('d', markerSymbol)
                .attr('transform', d => `translate(${xScale(getVariantXPosition(d))},${variantTrackY}) rotate(180)`)
                .attr('fill', d => variantColorScale(d.type || 'VAR'))
                .on('mouseover', function(event, d) {
                    tooltip.transition()
                        .duration(100)
                        .style('opacity', .95);
                    tooltip.html(formatVariantTooltip(d))
                        .style('left', (event.pageX + 10) + 'px')
                        .style('top', (event.pageY - 28) + 'px');
                })
                .on('mouseout', function() {
                    tooltip.transition()
                        .duration(500)
                        .style('opacity', 0);
                })
                .on('click', function(event, d) {
                    highlightVariantRow(d.variant_index);
                });
        }

        function highlightVariantRow(variantIndex) {
            document.querySelectorAll('.variant-row-selected').forEach(row => {
                row.classList.remove('variant-row-selected');
            });

            const row = document.getElementById(`variant-row-${variantIndex}`);
            if (!row) return;
            row.classList.add('variant-row-selected');
            row.scrollIntoView({behavior: 'smooth', block: 'nearest'});
        }

        function updateVariantPanel(mappedVariants, currentGenomeValue) {
            const countBadge = document.getElementById('variant-count-badge');
            const fileList = document.getElementById('variant-file-list');
            const tableDiv = document.getElementById('variant-table');
            if (!countBadge || !fileList || !tableDiv) return;

            countBadge.textContent = `${mappedVariants.length.toLocaleString()} variant${mappedVariants.length === 1 ? '' : 's'}`;

            if (!vcfFiles.length) {
                fileList.innerHTML = `<span class="chip"><i class="pi pi-filter"></i>Min depth >= ${Number(minVariantDepth).toLocaleString()}</span>`;
                tableDiv.innerHTML = '<div class="empty-state">Use <code>--vcf sample.vcf.gz</code> with a matching <code>sample.vcf.gz.tbi</code> index to show variants.</div>';
                return;
            }

            fileList.innerHTML = `<span class="chip"><i class="pi pi-filter"></i>Min depth >= ${Number(minVariantDepth).toLocaleString()}</span><span class="chip"><i class="pi pi-file"></i>VCF source${vcfFiles.length === 1 ? '' : 's'}</span>${vcfFiles.map(f => `<span class="chip"><strong>${escapeHtml(f.name)}</strong></span>`).join('')}`;

            if (!mappedVariants.length) {
                const genomeName = currentGenomeValue && currentGenomeValue.name ? currentGenomeValue.name : 'this genome';
                tableDiv.innerHTML = `<div class="empty-state">No variants with depth >= ${Number(minVariantDepth).toLocaleString()} overlap the displayed signature fragments for ${escapeHtml(genomeName)}.</div>`;
                return;
            }

            const maxRows = 250;
            const rows = mappedVariants.slice(0, maxRows).map(d => `
                <tr id="variant-row-${d.variant_index}">
                    <td>${escapeHtml(d.seq_name || '')}</td>
                    <td>${Number(d.pos).toLocaleString()}${d.end && d.end !== d.pos ? '-' + Number(d.end).toLocaleString() : ''}</td>
                    <td>${escapeHtml(d.type || 'VAR')}</td>
                    <td>${escapeHtml(d.ref || '')}</td>
                    <td>${escapeHtml(d.alt || '')}</td>
                    <td>${escapeHtml(d.qual || '.')}</td>
                    <td>${formatNullableNumber(d.depth)}</td>
                    <td>${formatNullablePercent(d.ref_dp)}</td>
                    <td>${formatNullablePercent(d.alt_dp)}</td>
                </tr>
            `).join('');

            const truncatedNotice = mappedVariants.length > maxRows
                ? `<p class="text-muted small mb-2">Showing the first ${maxRows.toLocaleString()} of ${mappedVariants.length.toLocaleString()} overlapping variants. All variants are drawn on the plot.</p>`
                : '';

            tableDiv.innerHTML = `
                ${truncatedNotice}
                <div class="table-responsive variant-table-wrap">
                    <table class="table table-sm table-striped variant-table">
                        <thead>
                            <tr>
                                <th>Seq</th>
                                <th>Position</th>
                                <th>Type</th>
                                <th>REF</th>
                                <th>ALT</th>
                                <th>QUAL</th>
                                <th>DEPTH</th>
                                <th>REF_DP %</th>
                                <th>ALT_DP %</th>
                            </tr>
                        </thead>
                        <tbody>${rows}</tbody>
                    </table>
                </div>
            `;
        }
        
        function networkSpeciesKey(node) {
            return String(node.speciesName || '').trim();
        }

        function buildNetworkSpeciesRows(nodes) {
            return [
                ...new Map(
                    nodes.map(node => {
                        const speciesName = networkSpeciesKey(node);
                        const species = speciesData.find(
                            item => String(item.name) === speciesName
                        );

                        if (!species) return null;

                        return [
                            species.taxid,
                            {
                                taxid: species.taxid,
                                name: speciesName,
                                readCount: species.readCount,
                                sigCov: species.sigCov,
                                sniScore: species.sniScore,
                            }
                        ];
                    }).filter(Boolean)
                ).values()
            ].sort((a, b) => (b.readCount ?? 0) - (a.readCount ?? 0));
        }

        function updateNetworkSpeciesTable(nodes, selectedId) {
            const container = document.getElementById('network-species-table');
            if (!container) return;

            if (!nodes.length) {
                container.innerHTML = '<div class="empty-state">No relevant species to display.</div>';
                return;
            }

            const rows = buildNetworkSpeciesRows(nodes).map(row => `
                <tr>
                    <td>${escapeHtml(row.name)}</td>
                    <td>${formatProfileValue(row.readCount, 'count')}</td>
                    <td>${formatProfileValue(row.sigCov, 'percent')}</td>
                    <td>${formatProfileValue(row.sniScore, 'score')}</td>
                </tr>
            `).join('');

            container.innerHTML = `
                <div class="table-responsive network-species-table-wrap">
                    <table class="table table-sm table-striped network-species-table">
                        <thead>
                            <tr>
                                <th>Species</th>
                                <th>Read Count</th>
                                <th>Signature Coverage</th>
                                <th>SNI Score</th>
                            </tr>
                        </thead>
                        <tbody>${rows}</tbody>
                    </table>
                </div>
            `;
        }

        const allNetworkSpecies = [...new Set(networkNodes.map(networkSpeciesKey).filter(Boolean))].sort();
        const elegantNetworkPalette = [
            '#52789c', '#69a297', '#8b78a6', '#c28b62', '#6f92b5',
            '#8ea66f', '#b87583', '#5f9ea0', '#9a7f6f', '#7f8db3',
            '#78a6a3', '#a68a69', '#8d78a8', '#6f9a7e', '#b47d6d',
            '#7189a8', '#9a9271', '#8a7f99', '#6f9b9a', '#a67f8d'
        ];
        const networkSpeciesColors = allNetworkSpecies.map((species, index) =>
            elegantNetworkPalette[index % elegantNetworkPalette.length]
        );
        const networkSpeciesColor = d3.scaleOrdinal()
            .domain(allNetworkSpecies)
            .range(networkSpeciesColors);

        function networkNodeColor(node) {
            const species = networkSpeciesKey(node);
            return species ? networkSpeciesColor(species) : '#a5adb8';
        }

        function updateNetworkGraph(currentGenomeValue) {
            const plot = document.getElementById('network-plot');
            const legend = d3.select('#network-legend').html('');
            const badge = document.getElementById('network-count-badge');
            if (!plot || !badge) return;

            if (networkSimulation) {
                networkSimulation.stop();
                networkSimulation = null;
            }
            d3.select(plot).html('');

            if (!currentGenomeValue || !networkNodes.length) {
                badge.textContent = '0 nodes · 0 edges';
                plot.innerHTML = '<div class="empty-state">No supplementary-alignment network data are available.</div>';
                updateNetworkSpeciesTable([]);
                return;
            }

            const selectedId = String(currentGenomeValue.taxid);
            const relevantLinks = networkEdges
                .filter(edge => String(edge.source) === selectedId || String(edge.target) === selectedId)
                .map(edge => ({ ...edge, source: String(edge.source), target: String(edge.target) }));

            // Always show the selected genome's node, even when it has no cross-TaxID edges.
            const relevantIds = new Set([selectedId]);
            relevantLinks.forEach(edge => {
                relevantIds.add(String(edge.source));
                relevantIds.add(String(edge.target));
            });

            const nodeByTaxid = new Map(networkNodes.map(node => [String(node.taxid), node]));
            const nodes = [...relevantIds].map(taxid => {
                const source = nodeByTaxid.get(taxid) || (taxid === selectedId
                    ? { taxid: selectedId, name: currentGenomeValue.name, speciesName: currentGenomeValue.parentName || '' }
                    : { taxid });
                return { ...source, id: taxid };
            });
            const nodeIds = new Set(nodes.map(node => node.id));
            const links = relevantLinks.filter(edge => nodeIds.has(edge.source) && nodeIds.has(edge.target));

            badge.textContent = `${nodes.length.toLocaleString()} node${nodes.length === 1 ? '' : 's'} · ${links.length.toLocaleString()} edge${links.length === 1 ? '' : 's'}`;
            updateNetworkSpeciesTable(nodes, selectedId);

            const speciesShown = [...new Set(nodes.map(networkSpeciesKey).filter(Boolean))].sort();
            speciesShown.forEach(species => {
                const item = legend.append('div').attr('class', 'legend-item');
                item.append('div').attr('class', 'legend-color').style('background-color', networkSpeciesColor(species));
                item.append('span').text(species);
            });
            if (nodes.some(node => !networkSpeciesKey(node))) {
                const item = legend.append('div').attr('class', 'legend-item');
                item.append('div').attr('class', 'legend-color').style('background-color', '#a5adb8');
                item.append('span').text('Species unavailable');
            }

            const outerWidth = Math.max(640, Math.floor(plot.getBoundingClientRect().width || 900));
            const outerHeight = 520;
            const graphSvg = d3.select(plot)
                .append('svg')
                .attr('viewBox', `0 0 ${outerWidth} ${outerHeight}`)
                .attr('preserveAspectRatio', 'xMidYMid meet')
                .attr('role', 'img')
                .attr('aria-label', `Supplementary alignment network for TaxID ${selectedId}`);

            const defs = graphSvg.append('defs');
            defs.append('marker')
                .attr('id', 'network-arrow')
                .attr('viewBox', '0 -4 8 8')
                .attr('refX', 7)
                .attr('refY', 0)
                .attr('markerWidth', 3)
                .attr('markerHeight', 3)
                .attr('orient', 'auto')
                .append('path')
                .attr('d', 'M0,-4L8,0L0,4')
                .attr('fill', '#788697');

            const maxNodeReads = d3.max(nodes, node => Number(node.readCount || 0)) || 1;
            const radiusScale = d3.scaleSqrt().domain([0, maxNodeReads]).range([8, 22]);
            const maxEdgeReads = d3.max(links, edge => Number(edge.readCount || 0)) || 1;
            const edgeWidthScale = d3.scaleSqrt().domain([1, maxEdgeReads]).range([1.2, 6]);
            const nodeRadius = node => radiusScale(Number(node.readCount || 0)) + (node.id === selectedId ? 3 : 0);

            const pairCounts = new Map();
            links.forEach(edge => {
                const pair = [edge.source, edge.target].sort().join('|');
                pairCounts.set(pair, (pairCounts.get(pair) || 0) + 1);
            });

            const link = graphSvg.append('g')
                .selectAll('path')
                .data(links)
                .enter()
                .append('path')
                .attr('class', 'network-link')
                .attr('stroke-width', edge => edgeWidthScale(Math.max(1, Number(edge.readCount || 0))))
                .attr('marker-end', 'url(#network-arrow)')
                .on('mouseover', function(event, edge) {
                    tooltip.transition().duration(100).style('opacity', .95);
                    tooltip.html(`
                        <strong>${escapeHtml(edge.source.id || edge.source)} → ${escapeHtml(edge.target.id || edge.target)}<br>
                        <strong>Reads:</strong> ${Number(edge.readCount || 0).toLocaleString()}<br>
                        <strong>Alignments:</strong> ${Number(edge.alignmentCount || 0).toLocaleString()}
                    `)
                    .style('left', (event.pageX + 10) + 'px')
                    .style('top', (event.pageY - 28) + 'px');
                })
                .on('mouseout', function() {
                    tooltip.transition().duration(300).style('opacity', 0);
                });

            const nodeGroup = graphSvg.append('g')
                .selectAll('g')
                .data(nodes)
                .enter()
                .append('g');

            nodeGroup.append('circle')
                .attr('class', node => `network-node${node.id === selectedId ? ' selected' : ''}${genomeData.some(genome => String(genome.taxid) === node.id) ? ' clickable' : ''}`)
                .attr('r', nodeRadius)
                .attr('fill', networkNodeColor)
                .on('mouseover', function(event, node) {
                    tooltip.transition().duration(100).style('opacity', .95);
                    tooltip.html(`
                        <strong>${escapeHtml(node.name || `TaxID ${node.id}`)}</strong><br>
                        <strong>Species:</strong> ${escapeHtml(node.speciesName || 'N/A')}<br>
                        <strong>TaxID:</strong> ${escapeHtml(node.id)}<br>
                        <strong>Reads:</strong> ${Number(node.readCount || 0).toLocaleString()}<br>
                        <strong>Alignments:</strong> ${Number(node.alignmentCount || 0).toLocaleString()}<br>
                        <strong>Primary:</strong> ${Number(node.primaryAlignmentCount || 0).toLocaleString()}<br>
                        <strong>Supplementary:</strong> ${Number(node.suppAlignmentCount || 0).toLocaleString()}
                    `)
                    .style('left', (event.pageX + 10) + 'px')
                    .style('top', (event.pageY - 28) + 'px');
                })
                .on('mouseout', function() {
                    tooltip.transition().duration(300).style('opacity', 0);
                })
                .on('click', function(event, node) {
                    if (node.id === selectedId) return;
                    const genome = genomeData.find(item => String(item.taxid) === node.id);
                    if (!genome) return;
                    selectedTaxid.value = node.id;
                    selectedGenome.value = node.id;
                    currentGenome.value = genome;
                    updateVisualization();
                });

            nodeGroup.append('text')
                .attr('class', node => `network-label${node.id === selectedId ? ' selected' : ''}`)
                .attr('text-anchor', 'middle')
                .attr('dy', node => nodeRadius(node) + 14)
                .text(node => {
                    const label = node.name || `TaxID ${node.id}`;
                    return label.length > 28 ? label.slice(0, 25) + '…' : label;
                });

            function dragstarted(event, node) {
                if (!event.active) networkSimulation.alphaTarget(0.3).restart();
                node.fx = node.x;
                node.fy = node.y;
            }

            function dragged(event, node) {
                node.fx = event.x;
                node.fy = event.y;
            }

            function dragended(event, node) {
                if (!event.active) networkSimulation.alphaTarget(0);
                node.fx = null;
                node.fy = null;
            }

            nodeGroup.call(d3.drag()
                .on('start', dragstarted)
                .on('drag', dragged)
                .on('end', dragended));

            networkSimulation = d3.forceSimulation(nodes)
                .force('link', d3.forceLink(links).id(node => node.id).distance(125).strength(0.7))
                .force('charge', d3.forceManyBody().strength(-420))
                .force('center', d3.forceCenter(outerWidth / 2, outerHeight / 2))
                .force('collision', d3.forceCollide().radius(node => nodeRadius(node) + 20))
                .on('tick', () => {
                    nodes.forEach(node => {
                        const r = nodeRadius(node) + 30;
                        node.x = Math.max(r, Math.min(outerWidth - r, node.x));
                        node.y = Math.max(r, Math.min(outerHeight - r, node.y));
                    });

                    link.attr('d', edge => {
                        const source = edge.source;
                        const target = edge.target;
                        const dx = target.x - source.x;
                        const dy = target.y - source.y;
                        const distance = Math.max(1, Math.sqrt(dx * dx + dy * dy));
                        const ux = dx / distance;
                        const uy = dy / distance;
                        const sourceR = nodeRadius(source) + 2;
                        const targetR = nodeRadius(target) + 6;
                        const sx = source.x + ux * sourceR;
                        const sy = source.y + uy * sourceR;
                        const tx = target.x - ux * targetR;
                        const ty = target.y - uy * targetR;
                        const pair = [source.id, target.id].sort().join('|');
                        if ((pairCounts.get(pair) || 0) > 1) {
                            const mx = (sx + tx) / 2 - uy * 24;
                            const my = (sy + ty) / 2 + ux * 24;
                            return `M${sx},${sy} Q${mx},${my} ${tx},${ty}`;
                        }
                        return `M${sx},${sy} L${tx},${ty}`;
                    });
                    nodeGroup.attr('transform', node => `translate(${node.x},${node.y})`);
                });
        }

        // Update the visualization for the selected genome
        function updateVisualization() {
            initVisualization();

            console.log("Updating visualization for genome:", currentGenome.value);
            
            // Access currentGenome through the global reactive reference
            const currentGenomeValue = currentGenome.value;
            if (!currentGenomeValue) return;
            
            const taxid = currentGenomeValue.taxid;
            updateNetworkGraph(currentGenomeValue);
            const genomeFragments = coverageData.filter(d => d.genome_taxid === taxid);
            
            if (genomeFragments.length === 0) {
                d3.select("#coverage-plot").html("<p>No coverage data available for this genome.</p>");
                updateVariantPanel([], currentGenomeValue);
                return;
            }
            
            // Sort fragments by sequence name and position
            genomeFragments.sort((a, b) => {
                if (a.seq_name === b.seq_name) {
                    return a.start_position - b.start_position;
                }
                return a.seq_name.localeCompare(b.seq_name);
            });
            
            // Set the x domain for the currently selected coordinate system
            xScale.domain(getFullXAxisDomain(currentGenomeValue));
            
            // Set y domain based on maximum mean depth with some padding
            const maxMeanDepth = d3.max(genomeFragments, d => d.meandepth);
            yScale.domain([0, maxMeanDepth * 1.2]); // Add 20% padding
            
            // Create percentage scale for top axis (0-100%)
            const xScalePercent = d3.scaleLinear()
                .domain([0, 100])
                .range([0, width]);
            
            // Update axes
            xAxis.call(d3.axisBottom(xScale));
            xAxisTop.call(d3.axisTop(xScalePercent).tickFormat(d => d + "%"));
            yAxis.call(d3.axisLeft(yScale));
            
            // Create color scale for different sequences
            const seqNames = [...new Set(genomeFragments.map(d => d.seq_name))];
            const colorScale = d3.scaleOrdinal(d3.schemeTableau10 || d3.schemeCategory10)
                .domain(seqNames);
            
            // Create legend
            const legend = d3.select("#legend").html("");
            seqNames.forEach(seq => {
                const item = legend.append("div")
                    .attr("class", "legend-item");
                
                item.append("div")
                    .attr("class", "legend-color")
                    .style("background-color", colorScale(seq));
                
                item.append("span")
                    .text(seq);
            });
            
            // Calculate positions - concatenate fragments without gaps
            let position = 0;
            genomeFragments.forEach(fragment => {
                fragment.x_start = position;
                const fragmentLength = fragment.end_position - fragment.start_position + 1;
                fragment.x_end = position + fragmentLength;
                position = fragment.x_end; // Concatenate signature fragments without gaps
            });
            
            // Draw fragments - uncovered portion in gray, covered portion in color
            // First draw the entire fragment in gray
            svg.selectAll(".fragment-background")
                .data(genomeFragments)
                .enter()
                .append("rect")
                .attr("class", "fragment-background")
                .attr("x", d => xScale(getFragmentXStart(d)))
                .attr("y", d => yScale(d.meandepth))
                .attr("width", getFragmentPixelWidth)
                .attr("height", d => height - yScale(d.meandepth))
                .attr("fill", "#e8eef6")
                .attr("rx", 5)
                .attr("ry", 5);
            
            // Add signature fragment indicators on top of histogram
            svg.selectAll(".fragment-indicator")
                .data(genomeFragments)
                .enter()
                .append("rect")
                .attr("class", "fragment-indicator")
                .attr("x", d => xScale(getFragmentXStart(d)))
                .attr("y", 0) // Position at the top below the axis
                .attr("width", getFragmentPixelWidth)
                .attr("height", 15)
                .attr("fill", d => colorScale(d.seq_name))
                .attr("opacity", d => d.coverage / 100) // Opacity based on coverage percentage
                .attr("rx", 4)
                .attr("ry", 4)
                .on("mouseover", function(event, d) {
                    tooltip.transition()
                        .duration(100)
                        .style("opacity", .95);
                    tooltip.html(`
                        <strong>Sequence:</strong> ${d.seq_name}<br>
                        <strong>Position:</strong> ${d.start_position}-${d.end_position}<br>
                        <strong>Coverage:</strong> ${d.coverage.toFixed(2)}%<br>
                        <strong>Read count:</strong> ${d.numreads.toFixed(0)}<br>
                        <strong>Mean Depth:</strong> ${d.meandepth.toFixed(2)}x<br>
                        <small class="text-muted">Click to view NCBI annotations</small>
                    `)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                })
                .on("mouseout", function() {
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                })
                .on("click", function(event, d) {
                    // Store current fragment data for NCBI query in Vue's reactive state
                    currentFragment.value = d;
                    
                    // Show the modal with loading indicator
                    document.getElementById('annotation-results').innerHTML = '';
                    document.getElementById('annotation-loading').classList.remove('d-none');
                    annotationModal.show();
                    
                    // Set modal title with fragment info
                    document.getElementById('annotationModalLabel').textContent = 
                        `NCBI Annotations: ${d.seq_name} (${d.start_position}-${d.end_position})`;
                    
                    // Query NCBI for annotations
                    queryNucleotide(d.seq_name, d.start_position, d.end_position);
                });
            
            // Then draw the covered portion in color
            svg.selectAll(".fragment")
                .data(genomeFragments)
                .enter()
                .append("rect")
                .attr("class", "fragment")
                .attr("x", d => xScale(getFragmentXStart(d)))
                .attr("y", d => yScale(d.meandepth))
                .attr("width", d => getFragmentPixelWidth(d) * (d.coverage / 100)) // Width based on coverage percentage
                .attr("height", d => height - yScale(d.meandepth))
                .attr("fill", d => colorScale(d.seq_name))
                .attr("rx", 4)
                .attr("ry", 3)
                .on("mouseover", function(event, d) {
                    tooltip.transition()
                        .duration(100)
                        .style("opacity", .95);
                    tooltip.html(`
                        <strong>Sequence:</strong> ${d.seq_name}<br>
                        <strong>Position:</strong> ${d.start_position}-${d.end_position}<br>
                        <strong>Coverage:</strong> ${d.coverage.toFixed(2)}%<br>
                        <strong>Read count:</strong> ${d.numreads.toFixed(0)}<br>
                        <strong>Mean Depth:</strong> ${d.meandepth.toFixed(2)}x<br>
                        <small class="text-muted">Click to view NCBI annotations</small>
                    `)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                })
                .on("mouseout", function() {
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                })
                .on("click", function(event, d) {
                    // Store current fragment data for NCBI query in Vue's reactive state
                    currentFragment.value = d;
                    
                    // Show the modal with loading indicator
                    document.getElementById('annotation-results').innerHTML = '';
                    document.getElementById('annotation-loading').classList.remove('d-none');
                    annotationModal.show();
                    
                    // Set modal title with fragment info
                    document.getElementById('annotationModalLabel').textContent = 
                        `NCBI Annotations: ${d.seq_name} (${d.start_position}-${d.end_position})`;
                    
                    // Query NCBI for annotations
                    queryNucleotide(d.seq_name, d.start_position, d.end_position);
                });

            
            // Draw mean depth line
            const overallMeanDepth = d3.mean(genomeFragments, d => d.meandepth);
            svg.append("line")
                .attr("x1", 0)
                .attr("y1", yScale(overallMeanDepth))
                .attr("x2", width)
                .attr("y2", yScale(overallMeanDepth))
                .attr("stroke", "#be3455")
                .attr("stroke-width", 2)
                .attr("stroke-opacity", 0.72)
                .attr("stroke-dasharray", "5,5")
                .style("cursor", "pointer")
                .on("mouseover", function(event) {
                    tooltip.transition()
                        .duration(200)
                        .style("opacity", .9);
                    tooltip.html(`<strong>Overall Mean Depth:</strong> ${overallMeanDepth.toFixed(2)}x`)
                        .style("left", (event.pageX + 10) + "px")
                        .style("top", (event.pageY - 28) + "px");
                })
                .on("mouseout", function() {
                    tooltip.transition()
                        .duration(500)
                        .style("opacity", 0);
                });

            const mappedVariants = getMappedVariants(genomeFragments);
            drawVariants(mappedVariants);
            updateVariantPanel(mappedVariants, currentGenomeValue);
        }
        
        // NCBI Nucleotide Annotation Query function
        async function queryNucleotide(accession, startPosition, stopPosition) {
            const resultsDiv = document.getElementById('annotation-results');
            const loadingDiv = document.getElementById('annotation-loading');
            
            try {
                // Check if accession is a generic sequence descriptor, not a real accession
                const genericTerms = ["chromosome", "scaffold", "contig", "node", "chr"];
                const accessionLower = accession.toLowerCase();
                
                for (const term of genericTerms) {
                    if (accessionLower.startsWith(term)) {
                        resultsDiv.innerHTML = `<div class="alert alert-info">
                            <i class="bi bi-info-circle me-2"></i>
                            Sequence "${accession}" appears to be a generic sequence identifier rather than an NCBI accession.
                            NCBI annotations are not available for this sequence.
                        </div>`;
                        loadingDiv.classList.add('d-none');
                        return;
                    }
                }
                
                // Set up NCBI API URLs
                const baseUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
                const encodedAccession = encodeURIComponent(accession);

                // Step 1: Search for the accession to get the UID in the 'nuccore' database
                const esearchUrl = `${baseUrl}esearch.fcgi?db=nuccore&term=${encodedAccession}&retmode=json`;
                console.log('ESearch URL:', esearchUrl);
                const esearchResponse = await fetch(esearchUrl);
                const esearchData = await esearchResponse.json();

                if (esearchData.esearchresult.idlist && esearchData.esearchresult.idlist.length > 0) {
                    const uid = esearchData.esearchresult.idlist[0];

                    // Step 2: Fetch the nucleotide record in GenBank format
                    const efetchUrl = `${baseUrl}efetch.fcgi?db=nuccore&id=${encodeURIComponent(uid)}&seq_start=${encodeURIComponent(startPosition)}&seq_stop=${encodeURIComponent(stopPosition)}&rettype=gb&retmode=text`;
                    console.log('EFetch URL:', efetchUrl);
                    const efetchResponse = await fetch(efetchUrl);
                    const gbText = await efetchResponse.text();
                    
                    // Step 3: Parse the GenBank format output
                    const lines = gbText.split('\n');
                    let sequenceName = '';
                    let organismName = '';
                    let definition = '';
                    let plasmid = 'No';
                    const annotations = [];
                    let currentFeature = null;
                    let inFeatures = false;
                    
                    // Estimate sequence length from the query range
                    const sequenceLength = stopPosition - startPosition + 1;
                    
                    // Parse the GenBank header information and features
                    for (let i = 0; i < lines.length; i++) {
                        const line = lines[i].trim();
                        
                        // Parse header information
                        if (line.startsWith('LOCUS')) {
                            // Parse LOCUS line for sequence name
                            const parts = line.split(/\s+/);
                            if (parts.length > 1) {
                                sequenceName = parts[1];
                            }
                        }
                        else if (line.startsWith('DEFINITION')) {
                            // Get sequence definition
                            definition = line.substring(10).trim();
                            // Continue to collect multi-line definitions
                            let j = i + 1;
                            while (j < lines.length && !lines[j].trim().startsWith('ACCESSION') && lines[j].trim().length > 0) {
                                definition += ' ' + lines[j].trim();
                                j++;
                            }

                            plasmid = definition.toLowerCase().includes('plasmid') ? 'Yes' : 'No';
                        }
                        else if (line.startsWith('ORGANISM')) {
                            // Get organism name
                            organismName = line.substring(10).trim();
                        } 
                        // Start of features section
                        else if (line === 'FEATURES             Location/Qualifiers') {
                            inFeatures = true;
                            continue;
                        } 
                        // End of features section
                        else if (line.startsWith('ORIGIN')) {
                            inFeatures = false;
                            break;
                        }
                        // Parse features
                        else if (inFeatures) {
                            // Feature lines start with 5 spaces then the feature type (like "gene", "CDS", etc.)
                            if (line.match(/^[a-zA-Z]/) && !line.includes('Location/Qualifiers')) {
                                // Save the previous feature if it exists
                                if (currentFeature) {
                                    annotations.push(currentFeature);
                                }
                                
                                // Get feature type (first word)
                                const featureParts = line.split(/\s+/);
                                
                                // Create new feature
                                currentFeature = {
                                    type: featureParts[0],
                                    location: featureParts[1],
                                    qualifiers: {}
                                };
                            }
                            // Qualifier lines start with 21 spaces then /qualifier=
                            else if (line.match(/^\//) && currentFeature) {
                                // Extract qualifier name and value
                                const qualParts = line.split('=');
                                
                                if (qualParts.length > 1) {
                                    // Handle qualifier with value like /gene="lacZ"
                                    const qualifierName = qualParts[0].replace('/', ''); // Remove leading '/'
                                    let qualifierValue = qualParts[1];

                                    // Handle quoted values
                                    if (qualifierValue.startsWith('"')) {
                                        qualifierValue = qualifierValue.substring(1);
                                        
                                        // Check if the value ends with quote or continues on next lines
                                        if (qualifierValue.endsWith('"')) {
                                            qualifierValue = qualifierValue.substring(0, qualifierValue.length - 1);
                                        } else {
                                            // Multi-line qualifier value
                                            let j = i + 1;
                                            while (j < lines.length) {
                                                const nextLine = lines[j].trim();
                                                if (nextLine.endsWith('"')) {
                                                    qualifierValue += ' ' + nextLine.substring(0, nextLine.length - 1);
                                                    i = j; // Skip these lines in the main loop
                                                    break;
                                                } else {
                                                    qualifierValue += ' ' + nextLine;
                                                }
                                                j++;
                                            }
                                        }
                                    }
                                    
                                    // Add to feature's qualifiers. Skip if the qualifier name is translation
                                    if (qualifierName !== 'translation') {
                                        currentFeature.qualifiers[qualifierName] = qualifierValue;
                                    }
                                }
                            }
                        }
                    }
                    
                    // Add the last feature if it exists
                    if (currentFeature) {
                        annotations.push(currentFeature);
                    }

                    // Generate result HTML - sequence details in clean format
                    let output = `<div class="card mb-3 shadow-sm">
                        <div class="card-header text-white" style="background-color: #a7aca2;">
                            Sequence Details
                        </div>
                        <div class="card-body">
                            <div class="row">
                                <div class="col-md-6">
                                    <p><strong>Accession:</strong> ${accession}</p>
                                    <p><strong>Name:</strong> ${sequenceName || 'N/A'}</p>
                                    <p><strong>Organism:</strong> ${organismName || 'N/A'}</p>
                                </div>
                                <div class="col-md-6">
                                    <p><strong>Length:</strong> ${sequenceLength.toLocaleString()} bp</p>
                                    <p><strong>Location:</strong> ${startPosition.toLocaleString()} - ${stopPosition.toLocaleString()}</p>
                                    <p><strong>Plasmid:</strong> ${plasmid || 'N/A'}</p>
                                </div>
                            </div>
                            <div class="row mt-2">
                                <div class="col-12">
                                    <p><strong>Definition:</strong> ${definition || 'N/A'}</p>
                                </div>
                            </div>
                        </div>
                    </div>`;

                    if (annotations.length > 0) {
                        // Group annotations by type for better organization
                        const groupedAnnotations = {};
                        annotations.forEach(anno => {
                            if (!groupedAnnotations[anno.type]) {
                                groupedAnnotations[anno.type] = [];
                            }
                            groupedAnnotations[anno.type].push(anno);
                        });
                        
                        output += `<div class="d-flex justify-content-between align-items-center mb-3">
                            <h5 class="mb-0">Annotations (${annotations.length})</h5>
                            <span class="badge" style="background-color: #b5a397; color: white">${startPosition.toLocaleString()}-${stopPosition.toLocaleString()} bp</span>
                        </div>`;
                        
                        // Create tabs for feature types
                        const featureTypes = Object.keys(groupedAnnotations);
                        
                        if (featureTypes.length > 1) {
                            // Add tab navigation
                            output += `<ul class="nav nav-tabs mb-3" id="annotationTabs" role="tablist">`;
                            featureTypes.forEach((type, index) => {
                                const isActive = index === 0 ? 'active' : '';
                                output += `
                                <li class="nav-item" role="presentation">
                                    <button class="nav-link ${isActive}" id="tab-${type}" data-bs-toggle="tab" 
                                    data-bs-target="#content-${type}" type="button" role="tab" aria-selected="${index === 0}">
                                        ${type} (${groupedAnnotations[type].length})
                                    </button>
                                </li>`;
                            });
                            output += `</ul>`;
                            
                            // Add tab content
                            output += `<div class="tab-content" id="annotationTabContent">`;
                            featureTypes.forEach((type, index) => {
                                const isActive = index === 0 ? 'show active' : '';
                                output += `
                                <div class="tab-pane fade ${isActive}" id="content-${type}" role="tabpanel" aria-labelledby="tab-${type}">`;
                                
                                // Feature cards in each tab
                                groupedAnnotations[type].forEach(anno => {
                                    output += createAnnotationCard(anno);
                                });
                                
                                output += `</div>`;
                            });
                            output += `</div>`;
                        } else {
                            // If only one feature type, no need for tabs
                            const type = featureTypes[0];
                            const annos = groupedAnnotations[type];
                            
                            annos.forEach(anno => {
                                output += createAnnotationCard(anno);
                            });
                        }
                    } else {
                        output += `<div class="alert alert-info mt-3">
                            <i class="bi bi-info-circle me-2"></i>
                            No specific annotations found overlapping the range ${startPosition.toLocaleString()}-${stopPosition.toLocaleString()}.
                        </div>`;
                    }
                    
                    // Helper function to create an annotation card
                    function createAnnotationCard(anno) {
                        let result = `
                        <div class="card mb-3 shadow-sm annotation-card">
                            <div class="card-body">
                                <div class="d-flex justify-content-between align-items-start">
                                    <h5 class="text-dark">${anno.type}</h5>
                                    <span class="badge bg-light text-dark">${anno.location || anno.start + '..' + anno.end}</span>
                                </div>`;
                        
                        // Display important qualifiers at the top
                        const priorityQualifiers = ['gene', 'product', 'note', 'function', 'protein_id'];
                        
                        if (Object.keys(anno.qualifiers).length > 0) {
                            // First show priority qualifiers
                            result += `<div class="annotation-qualifiers mt-2">`;
                            
                            // Show priority qualifiers first
                            priorityQualifiers.forEach(key => {
                                if (anno.qualifiers[key]) {
                                    const value = anno.qualifiers[key];
                                    if (value === true) {
                                        result += `<div class="mb-1"><span class="badge bg-secondary">${key}</span></div>`;
                                    } else {
                                        result += `<div class="mb-1"><strong>${key}:</strong> ${value}</div>`;
                                    }
                                }
                            });
                            
                            // Then show all other qualifiers
                            for (const [key, value] of Object.entries(anno.qualifiers)) {
                                // Skip already shown priority qualifiers
                                if (priorityQualifiers.includes(key)) continue;
                                
                                if (value === true) {
                                    result += `<div class="mb-1"><span class="badge bg-secondary">${key}</span></div>`;
                                } else {
                                    result += `<div class="mb-1"><strong>${key}:</strong> ${value}</div>`;
                                }
                            }
                            
                            result += `</div>`;
                        }
                        
                        result += `</div></div>`;
                        return result;
                    }

                    // Update results
                    resultsDiv.innerHTML = output;
                } else {
                    resultsDiv.innerHTML = `<div class="alert alert-danger">
                        Nucleotide sequence with accession "${accession}" not found.
                    </div>`;
                }
            } catch (error) {
                console.error('Error fetching data:', error);
                resultsDiv.innerHTML = `<div class="alert alert-danger">
                    An error occurred: ${error.message}. Please check the accession number and your internet connection.
                </div>`;
            } finally {
                loadingDiv.classList.add('d-none');
            }
        }
        
        // Update genome information panel
        function updateGenomeInfo() {
            if (!currentGenome) return;
            
            const infoDiv = document.getElementById('genome-info');
            infoDiv.innerHTML = `
                <div class="row">
                    <div class="col-md-6">
                        <p><strong>Name:</strong> ${currentGenome.name}</p>
                        <p><strong>Taxid:</strong> ${currentGenome.taxid}</p>
                        <p><strong>Domain:</strong> ${currentGenome.superkingdom}</p>
                    </div>
                    <div class="col-md-6">
                        <p><strong>Number of Sequences:</strong> ${currentGenome.numOfSeq}</p>
                        <p><strong>Total Length:</strong> ${currentGenome.totalLength.toLocaleString()}</p>
                        <p><strong>Genome Size:</strong> ${currentGenome.genomeSize.toLocaleString()}</p>
                    </div>
                </div>
            `;
        }
        
        // Handle brush zoom
        function brushended(event) {
            if (!event.selection) return;
            const [x0, x1] = event.selection.map(xScale.invert);
            xScale.domain([x0, x1]);
            
            // Remove brush
            svg.select(".brush").call(brush.move, null);
            
            // Update chart with new domain
            updateChart();
        }
        
        // Reset zoom
        function resetZoom() {
            const currentGenomeValue = currentGenome.value;
            if (!currentGenomeValue) return;
            xScale.domain(getFullXAxisDomain(currentGenomeValue));
            svg.select(".brush").call(brush.move, null);
            updateChart();
        }
        
        // Update chart after zoom
        function updateChart() {
            // Update x-axis
            svg.select(".x-axis")
                .interrupt()
                .transition()
                .duration(chartTransitionDuration)
                .ease(d3.easeCubicInOut)
                .call(d3.axisBottom(xScale));
            
            // Update top x-axis (percentage scale needs to be recalculated based on current zoom)
            const currentDomain = xScale.domain();
            const fullDomain = getFullXAxisDomain(currentGenome.value);
            const axisLength = fullDomain[1] - fullDomain[0];
            const startPercent = ((currentDomain[0] - fullDomain[0]) / axisLength) * 100;
            const endPercent = ((currentDomain[1] - fullDomain[0]) / axisLength) * 100;
            
            const xScalePercentZoomed = d3.scaleLinear()
                .domain([startPercent, endPercent])
                .range([0, width]);
            
            svg.select(".x-axis-top")
                .interrupt()
                .transition()
                .duration(chartTransitionDuration)
                .ease(d3.easeCubicInOut)
                .call(d3.axisTop(xScalePercentZoomed).tickFormat(d => d.toFixed(1) + "%"));

            svg.select(".x-axis-label")
                .interrupt()
                .attr("opacity", 0)
                .text(getBottomXAxisLabel())
                .transition()
                .duration(chartTransitionDuration / 2)
                .attr("opacity", 1);

            svg.select(".x-axis-top-label")
                .interrupt()
                .attr("opacity", 0)
                .text(getTopXAxisLabel())
                .transition()
                .duration(chartTransitionDuration / 2)
                .attr("opacity", 1);
            
            // Update fragments - all components including indicators
            svg.selectAll(".fragment-background")
                .interrupt()
                .transition()
                .duration(chartTransitionDuration)
                .ease(d3.easeCubicInOut)
                .attr("x", d => xScale(getFragmentXStart(d)))
                .attr("width", getFragmentPixelWidth);
            
            svg.selectAll(".fragment-indicator")
                .interrupt()
                .transition()
                .duration(chartTransitionDuration)
                .ease(d3.easeCubicInOut)
                .attr("x", d => xScale(getFragmentXStart(d)))
                .attr("width", getFragmentPixelWidth);
                
            svg.selectAll(".fragment")
                .interrupt()
                .transition()
                .duration(chartTransitionDuration)
                .ease(d3.easeCubicInOut)
                .attr("x", d => xScale(getFragmentXStart(d)))
                .attr("width", d => getFragmentPixelWidth(d) * (d.coverage / 100));

            const [domainStart, domainEnd] = xScale.domain();
            svg.selectAll(".variant-marker")
                .interrupt()
                .style("display", d => {
                    const position = getVariantXPosition(d);
                    return position >= domainStart && position <= domainEnd ? null : "none";
                })
                .transition()
                .duration(chartTransitionDuration)
                .ease(d3.easeCubicInOut)
                .attr("transform", d => `translate(${xScale(getVariantXPosition(d))},24) rotate(180)`);
        }
    }).catch(function(error) {
        console.error('Unable to load embedded report resources:', error);
        const appElement = document.getElementById('app');
        appElement.innerHTML = '<div class="empty-state">Unable to load this report. Please open it in a current version of your browser.</div>';
        appElement.style.visibility = 'visible';
        const pageLoader = document.getElementById('page-loader');
        if (pageLoader) pageLoader.remove();
    });
    </script>
</body>
</html>
"""

def parse_coverage_file(coverage_file):
    """Parse the coverage file and extract relevant information using pandas."""
    # Read the coverage file with pandas
    # Define column names based on the file format
    col_names = ['rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq']
    df = pd.read_csv(coverage_file, sep='\t', comment='#', names=col_names, header=None)
    
    # Parse the rname field which has format: [seq_name]|[start_position]|[stop_position]|[genom_taxid]
    df[['seq_name', 'start_position', 'end_position', 'genome_taxid']] = df['rname'].str.rstrip('|').str.split('|', expand=True)
    
    # Convert data types for numerical columns
    numeric_cols = ['startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq']
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric)
    df['start_position'] = pd.to_numeric(df['start_position'], errors='coerce').astype('Int64')
    df['end_position'] = pd.to_numeric(df['end_position'], errors='coerce').astype('Int64')
    df['genome_taxid'] = df['genome_taxid'].astype(str)
    
    # Handle meanmapq - set to 0 if it's not a valid number
    df['meanmapq'] = pd.to_numeric(df['meanmapq'], errors='coerce').fillna(0)
    
    # Remove rows with missing data in essential columns
    df = df.dropna(subset=['seq_name', 'start_position', 'end_position', 'genome_taxid'])
    
    # Convert DataFrame to list of dictionaries for compatibility with the existing code
    coverage_data = df.to_dict('records')
    
    return coverage_data, df['genome_taxid'].unique().tolist()

PROFILE_COUNT_FIELDS = (
    'READ_COUNT', 'TOTAL_BP_MAPPED', 'COVERED_SIG_LEN', 'AOI_READ_COUNT',
    'TOTAL_READ_LEN', 'TOTAL_BP_MISMATCH', 'TOTAL_BP_INDEL', 'MAPPED_SIG_LEN',
    'TOTAL_SIG_LEN', 'GENOME_COUNT', 'GENOME_SIZE',
)
PROFILE_FLOAT_FIELDS = (
    'SNI_SCORE', 'SIG_COV', 'DEPTH', 'REL_ABUNDANCE_GC', 'REL_ABUNDANCE',
    'ALN_IDENTITY', 'CONSENSUS_SEQ_IDENTITY', 'SIG_COV_RAW', 'COVERED_SIG_DEPTH',
    'COVERED_MAPPED_SIG_COV', 'ZSCORE', 'GENOMIC_CONTENT_EST', 'ABUNDANCE',
    'REL_ABUNDANCE_DEPTH',
)
PROFILE_TEXT_FIELDS = (
    'LEVEL', 'NAME', 'TAXID', 'PARENT_NAME', 'PARENT_TAXID',
    'SNI_CI95_LH', 'SIG_LEVEL', 'NOTE',
)


def _profile_number(value, integer=False):
    """Keep unavailable metrics distinct from zero and safe to serialize as JSON."""
    if pd.isna(value) or not math.isfinite(value):
        return None
    return int(value) if integer else float(value)


def parse_full_file(full_file, uniq_taxid_list):
    """Parse taxonomy full.tsv data keyed by TAXID and return genome data."""
    if not full_file:
        return []

    df = pd.read_csv(
        full_file,
        sep='\t',
        comment='#',
        low_memory=False,
        dtype={field: 'string' for field in PROFILE_TEXT_FIELDS},
    )

    required_cols = {'TAXID', 'PARENT_NAME', 'SIG_COV', 'SIG_LEVEL', 'NAME', 'TOTAL_SIG_LEN', 'GENOME_SIZE'}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"Full taxonomy file is missing required columns: {', '.join(sorted(missing_cols))}")

    uniq_taxids = {str(t) for t in uniq_taxid_list}
    df['TAXID'] = df['TAXID'].astype(str)
    df_species = df[df['LEVEL'] == 'species'].copy()
    df = df[df['TAXID'].isin(uniq_taxids)]
    for field in PROFILE_COUNT_FIELDS + PROFILE_FLOAT_FIELDS:
        if field in df.columns:
            df[field] = pd.to_numeric(df[field], errors='coerce')

    # Convert to dictionary format for the JavaScript
    genome_data = []
    for _, row in df.iterrows():
        taxid = str(row['TAXID'])
        parent_name = '' if pd.isna(row['PARENT_NAME']) else str(row['PARENT_NAME'])
        profile = {
            field: '' if pd.isna(row[field]) else str(row[field])
            for field in PROFILE_TEXT_FIELDS if field in df.columns
        }
        for fields, integer in ((PROFILE_COUNT_FIELDS, True), (PROFILE_FLOAT_FIELDS, False)):
            profile.update({
                field: _profile_number(row[field], integer)
                for field in fields if field in df.columns
            })
        sig_cov = profile['SIG_COV']
        total_sig_len = profile['TOTAL_SIG_LEN'] or 0
        genome_size = profile['GENOME_SIZE'] or 0
        
        genome_data.append({
            'db_level': profile['SIG_LEVEL'],
            'name': profile['NAME'],
            'parentName': parent_name,
            'sigCov': sig_cov,
            'taxid': taxid,
            'superkingdom': row.get('SUPERKINGDOM', ''),
            'numOfSeq': int(row.get('NUM_FRAG', 0)) if 'NUM_FRAG' in df.columns and not pd.isna(row.get('NUM_FRAG')) else 0,
            'max': int(row.get('LONGEST_SIG_LEN', 0)) if 'LONGEST_SIG_LEN' in df.columns and not pd.isna(row.get('LONGEST_SIG_LEN')) else 0,
            'min': int(row.get('SHORTEST_SIG_LEN', 0)) if 'SHORTEST_SIG_LEN' in df.columns and not pd.isna(row.get('SHORTEST_SIG_LEN')) else 0,
            'totalLength': total_sig_len,
            'genomeSize': genome_size,
            'profile': profile,
        })

    species_data = []
    for _, row in df_species.iterrows():
        taxid = str(row['TAXID'])
        parent_name = '' if pd.isna(row['PARENT_NAME']) else str(row['PARENT_NAME'])
        profile = {
            field: '' if pd.isna(row[field]) else str(row[field])
            for field in PROFILE_TEXT_FIELDS if field in df_species.columns
        }
        for fields, integer in ((PROFILE_COUNT_FIELDS, True), (PROFILE_FLOAT_FIELDS, False)):
            profile.update({
                field: _profile_number(row[field], integer)
                for field in fields if field in df_species.columns
            })
        sig_cov = profile['SIG_COV']
        read_count = profile['READ_COUNT']
        sni_score = profile['SNI_SCORE']
        species_data.append({
            'db_level': profile['SIG_LEVEL'],
            'name': profile['NAME'],
            'parentName': parent_name,
            'sigCov': sig_cov,
            'taxid': taxid,
            'superkingdom': row.get('SUPERKINGDOM', ''),
            'readCount': read_count,
            'sniScore': sni_score
        })

    return genome_data, species_data


def _flatten_vcf_args(vcf_args):
    """Flatten --vcf arguments from argparse into a simple list."""
    if not vcf_args:
        return []

    flattened = []
    for item in vcf_args:
        if isinstance(item, (list, tuple)):
            flattened.extend(item)
        else:
            flattened.append(item)
    return flattened


def _find_tabix_index(vcf_file):
    """Return the expected .tbi path for a compressed VCF, or None if missing."""
    candidates = [f"{vcf_file}.tbi"]
    if vcf_file.endswith('.gz'):
        candidates.append(vcf_file[:-3] + '.tbi')

    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    return None


def _validate_vcf_inputs(vcf_files):
    for vcf_file in vcf_files:
        if not os.path.exists(vcf_file):
            raise FileNotFoundError(f"VCF file not found: {vcf_file}")
        if not vcf_file.endswith(('.vcf.gz', '.vcf.bgz', '.gz', '.bgz')):
            raise ValueError(
                f"VCF must be bgzip/gzip-compressed and tabix-indexed: {vcf_file}"
            )
        if _find_tabix_index(vcf_file) is None:
            raise FileNotFoundError(
                f"Missing tabix index for {vcf_file}. Expected {vcf_file}.tbi"
            )


def _truncate_text(value, max_length=500):
    if value is None:
        return ''
    value = str(value)
    if len(value) <= max_length:
        return value
    return value[: max_length - 3] + '...'


def _parse_vcf_info(info_value):
    if not info_value or info_value == '.':
        return {}

    info = {}
    for field in info_value.split(';'):
        if not field:
            continue
        if '=' in field:
            key, value = field.split('=', 1)
            info[key] = value
        else:
            info[field] = True
    return info


def _parse_int(value):
    if value is None or value is True or value == '':
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _parse_dp4(dp4_value):
    if not dp4_value or dp4_value is True:
        return None, None

    parts = str(dp4_value).split(',')
    if len(parts) < 4:
        return None, None

    values = [_parse_int(value) for value in parts[:4]]
    if any(value is None for value in values):
        return None, None

    ref_fwd, ref_rev, alt_fwd, alt_rev = values
    return ref_fwd + ref_rev, alt_fwd + alt_rev


def _depth_percent(count, depth, fallback_total):
    if count is None:
        return None

    denominator = depth if depth and depth > 0 else fallback_total
    if not denominator or denominator <= 0:
        return None

    return round((count / denominator) * 100, 1)


def _infer_variant_type(ref, alt, info):
    svtype = info.get('SVTYPE')
    if svtype and svtype is not True:
        return str(svtype)

    alt_alleles = [allele for allele in str(alt).split(',') if allele and allele != '.']
    if not alt_alleles:
        return 'REF'

    if any(allele.startswith('<') and allele.endswith('>') for allele in alt_alleles):
        return 'SV'
    if len(ref) == 1 and all(len(allele) == 1 for allele in alt_alleles):
        return 'SNV'
    if all(len(allele) == len(ref) for allele in alt_alleles):
        return 'MNV'
    if any(len(allele) > len(ref) for allele in alt_alleles) and any(len(allele) < len(ref) for allele in alt_alleles):
        return 'INDEL'
    if any(len(allele) > len(ref) for allele in alt_alleles):
        return 'INS'
    if any(len(allele) < len(ref) for allele in alt_alleles):
        return 'DEL'
    return 'COMPLEX'


def _merge_intervals(intervals):
    if not intervals:
        return []

    intervals = sorted((int(start), int(end)) if int(start) <= int(end) else (int(end), int(start)) for start, end in intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(interval) for interval in merged]


def _parse_signature_contig(contig):
    parts = str(contig).rstrip('|').split('|')
    if len(parts) < 4:
        return None

    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError:
        return None

    return {
        'seq_name': parts[0],
        'start': start,
        'end': end,
        'taxid': parts[3],
    }


def _build_interval_lookup(coverage_data):
    intervals_by_seq = defaultdict(list)
    fetch_intervals_by_contig = defaultdict(list)
    for row in coverage_data:
        seq_name = str(row.get('seq_name', ''))
        if not seq_name:
            continue
        start = int(row['start_position'])
        end = int(row['end_position'])
        intervals_by_seq[seq_name].append((start, end))

        contig = str(row.get('rname', ''))
        if contig:
            fetch_intervals_by_contig[contig].append((1, abs(end - start) + 1))
        fetch_intervals_by_contig[seq_name].append((start, end))

    merged_by_seq = {seq: _merge_intervals(intervals) for seq, intervals in intervals_by_seq.items()}
    starts_by_seq = {seq: [start for start, _ in intervals] for seq, intervals in merged_by_seq.items()}
    fetch_by_contig = {contig: _merge_intervals(intervals) for contig, intervals in fetch_intervals_by_contig.items()}
    return merged_by_seq, starts_by_seq, fetch_by_contig


def _overlaps_any_interval(seq_name, start, end, intervals_by_seq, starts_by_seq):
    intervals = intervals_by_seq.get(seq_name)
    if not intervals:
        return False

    starts = starts_by_seq[seq_name]
    candidate_index = bisect.bisect_right(starts, end) - 1
    return candidate_index >= 0 and intervals[candidate_index][1] >= start


def _variant_from_vcf_line(line, source_index, sample_names=None):
    parts = line.rstrip('\n').split('\t')
    if len(parts) < 8:
        return None

    chrom, pos, variant_id, ref, alt, qual, filter_value, info_value = parts[:8]
    try:
        local_pos = int(pos)
    except ValueError:
        return None

    info = _parse_vcf_info(info_value)
    end = info.get('END')
    try:
        local_end = int(end) if end and end is not True else local_pos + max(len(ref), 1) - 1
    except ValueError:
        local_end = local_pos + max(len(ref), 1) - 1

    signature_contig = _parse_signature_contig(chrom)
    if signature_contig:
        seq_name = signature_contig['seq_name']
        pos = signature_contig['start'] + local_pos - 1
        end = signature_contig['start'] + local_end - 1
    else:
        seq_name = chrom
        pos = local_pos
        end = local_end

    variant_type = _infer_variant_type(ref, alt, info)
    depth = _parse_int(info.get('DP'))
    ref_dp_count, alt_dp_count = _parse_dp4(info.get('DP4'))
    dp4_total = (ref_dp_count or 0) + (alt_dp_count or 0)
    ref_dp = _depth_percent(ref_dp_count, depth, dp4_total)
    alt_dp = _depth_percent(alt_dp_count, depth, dp4_total)

    return {
        'source_index': source_index,
        'seq_name': seq_name,
        'pos': pos,
        'end': end,
        'ref': _truncate_text(ref, 120),
        'alt': _truncate_text(alt, 120),
        'qual': '' if qual == '.' else qual,
        'type': variant_type,
        'depth': depth,
        'ref_dp': ref_dp,
        'alt_dp': alt_dp,
    }


def _read_vcf_sample_names(vcf_file):
    """Read VCF sample names from the #CHROM header line."""
    with gzip.open(vcf_file, 'rt', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('#CHROM'):
                return line.rstrip('\n').split('\t')[9:]
    return []


def _read_vcf_header_and_lines(vcf_file):
    """Yield VCF header data and variant lines using the standard library gzip module."""
    sample_names = []
    with gzip.open(vcf_file, 'rt', encoding='utf-8', errors='replace') as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header_parts = line.rstrip('\n').split('\t')
                sample_names = header_parts[9:]
                continue
            if line.startswith('#'):
                continue
            yield sample_names, line


def _read_vcf_lines_with_tabix(vcf_file, fetch_intervals_by_contig):
    """Yield VCF lines through the tabix index when pysam is available."""
    try:
        import pysam
    except ImportError:
        return None

    sample_names = _read_vcf_sample_names(vcf_file)

    def iterator():
        tabix_file = pysam.TabixFile(vcf_file)
        try:
            for contig, intervals in fetch_intervals_by_contig.items():
                for start, end in intervals:
                    try:
                        records = tabix_file.fetch(contig, max(0, int(start) - 1), int(end))
                    except ValueError:
                        continue
                    for line in records:
                        yield sample_names, line
        finally:
            tabix_file.close()

    return iterator()


def _iter_relevant_vcf_lines(vcf_file, fetch_intervals_by_contig):
    tabix_iter = _read_vcf_lines_with_tabix(vcf_file, fetch_intervals_by_contig)
    if tabix_iter is not None:
        try:
            yield from tabix_iter
            return
        except Exception as exc:
            print(
                f"Warning: could not read {vcf_file} through its tabix index ({exc}); "
                "falling back to a sequential gzip scan.",
                file=sys.stderr,
            )
    else:
        print(
            "Warning: pysam is not installed; scanning compressed VCF sequentially. "
            "Install pysam to use the .tbi index for faster VCF parsing.",
            file=sys.stderr,
        )

    yield from _read_vcf_header_and_lines(vcf_file)


def parse_vcf_files(vcf_files, coverage_data, min_depth=5):
    """Parse compressed, indexed VCF files and keep variants overlapping signature fragments."""
    vcf_files = _flatten_vcf_args(vcf_files)
    if not vcf_files:
        return [], []

    _validate_vcf_inputs(vcf_files)
    intervals_by_seq, starts_by_seq, fetch_intervals_by_contig = _build_interval_lookup(coverage_data)

    variant_data = []
    vcf_file_data = []

    for source_index, vcf_file in enumerate(vcf_files):
        source_name = os.path.basename(vcf_file)

        seen_variant_lines = set()
        for sample_names, line in _iter_relevant_vcf_lines(vcf_file, fetch_intervals_by_contig):
            if line in seen_variant_lines:
                continue
            seen_variant_lines.add(line)
            variant = _variant_from_vcf_line(line, source_index, sample_names)
            if variant is None:
                continue
            if len(vcf_files) == 1:
                variant.pop('source_index', None)
            if min_depth is not None and min_depth > 0:
                depth = variant.get('depth')
                if depth is None or depth < min_depth:
                    continue
            if not _overlaps_any_interval(
                variant['seq_name'], int(variant['pos']), int(variant['end']), intervals_by_seq, starts_by_seq
            ):
                continue
            variant_data.append(variant)

        vcf_file_data.append({'name': source_name})

    return variant_data, vcf_file_data


def parse_taxid_network(node_file, edge_file, full_file=None):
    """Parse TAXID node/edge TSVs and annotate nodes with NAME/PARENT_NAME."""
    if not node_file or not edge_file:
        return [], []

    node_file = Path(node_file)
    edge_file = Path(edge_file)
    if not node_file.exists() or not edge_file.exists():
        return [], []

    node_df = pd.read_csv(node_file, sep='\t', dtype={'TAXID': 'string'})
    edge_df = pd.read_csv(
        edge_file, sep='\t',
        dtype={'SOURCE_TAXID': 'string', 'TARGET_TAXID': 'string'},
    )

    metadata = {}
    if full_file and Path(full_file).exists():
        profile_df = pd.read_csv(
            full_file, sep='\t', comment='#', low_memory=False,
            usecols=lambda column: column in {'TAXID', 'NAME', 'PARENT_NAME'},
            dtype='string',
        )
        if 'TAXID' in profile_df.columns:
            profile_df['TAXID'] = profile_df['TAXID'].astype(str)
            for _, row in profile_df.drop_duplicates('TAXID').iterrows():
                taxid = str(row['TAXID'])
                metadata[taxid] = {
                    'name': '' if pd.isna(row.get('NAME')) else str(row.get('NAME')),
                    'speciesName': '' if pd.isna(row.get('PARENT_NAME')) else str(row.get('PARENT_NAME')),
                }

    numeric_node_cols = [
        'READ_COUNT', 'ALIGNMENT_COUNT', 'PRIMARY_ALIGNMENT_COUNT', 'SUPP_ALIGNMENT_COUNT'
    ]
    for column in numeric_node_cols:
        if column in node_df.columns:
            node_df[column] = pd.to_numeric(node_df[column], errors='coerce').fillna(0).astype(int)

    nodes = []
    for _, row in node_df.iterrows():
        taxid = str(row['TAXID'])
        meta = metadata.get(taxid, {})
        nodes.append({
            'taxid': taxid,
            'name': meta.get('name', ''),
            'speciesName': meta.get('speciesName', ''),
            'readCount': int(row.get('READ_COUNT', 0)),
            'alignmentCount': int(row.get('ALIGNMENT_COUNT', 0)),
            'primaryAlignmentCount': int(row.get('PRIMARY_ALIGNMENT_COUNT', 0)),
            'suppAlignmentCount': int(row.get('SUPP_ALIGNMENT_COUNT', 0)),
        })

    for column in ('READ_COUNT', 'ALIGNMENT_COUNT'):
        if column in edge_df.columns:
            edge_df[column] = pd.to_numeric(edge_df[column], errors='coerce').fillna(0).astype(int)

    edges = []
    for _, row in edge_df.iterrows():
        edges.append({
            'source': str(row['SOURCE_TAXID']),
            'target': str(row['TARGET_TAXID']),
            'readCount': int(row.get('READ_COUNT', 0)),
            'alignmentCount': int(row.get('ALIGNMENT_COUNT', 0)),
        })

    return nodes, edges


def _compact_coverage_for_html(coverage_data):
    compact = []
    for row in coverage_data:
        compact.append({
            'seq_name': str(row['seq_name']),
            'start_position': int(row['start_position']),
            'end_position': int(row['end_position']),
            'genome_taxid': str(row['genome_taxid']),
            'numreads': int(row['numreads']),
            'coverage': round(float(row['coverage']), 4),
            'meandepth': round(float(row['meandepth']), 4),
        })
    return compact


def _compact_json(data):
    # Report names and notes are text, including when they contain HTML delimiters.
    return json.dumps(data, separators=(',', ':')).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')


def generate_html(
    coverage_data,
    genome_data,
    species_data,
    output_file,
    node_file,
    edge_file,
    full_file=None,
    variant_data=None,
    vcf_file_data=None,
    min_variant_depth=5,
):
    """Generate the HTML file with the visualization."""
    variant_data = variant_data or []
    vcf_file_data = vcf_file_data or []
    network_nodes, network_edges = parse_taxid_network(node_file, edge_file, full_file)

    # Convert data to JSON for JavaScript
    coverage_json = _compact_json(_compact_coverage_for_html(coverage_data))
    genome_json = _compact_json(genome_data)
    species_json = _compact_json(species_data)
    variant_json = _compact_json(variant_data)
    vcf_file_json = _compact_json(vcf_file_data)
    network_node_json = _compact_json(network_nodes)
    network_edge_json = _compact_json(network_edges)

    html_content = _embed_browser_resources(HTML_TEMPLATE)

    # Replace placeholders in the HTML template
    html_content = html_content.replace('GENOME_DATA_PLACEHOLDER', genome_json)
    html_content = html_content.replace('SPECIES_DATA_PLACEHOLDER', species_json)
    html_content = html_content.replace('COVERAGE_DATA_PLACEHOLDER', coverage_json)
    html_content = html_content.replace('VARIANT_DATA_PLACEHOLDER', variant_json)
    html_content = html_content.replace('VCF_FILE_DATA_PLACEHOLDER', vcf_file_json)
    html_content = html_content.replace('NETWORK_NODE_DATA_PLACEHOLDER', network_node_json)
    html_content = html_content.replace('NETWORK_EDGE_DATA_PLACEHOLDER', network_edge_json)
    html_content = html_content.replace('MIN_VARIANT_DEPTH_PLACEHOLDER', str(min_variant_depth))
    
    if minify_html is not None:
        html_content = minify_html.minify(
            html_content,
            keep_comments=False,
            minify_css=True,
            minify_js=True,
            remove_processing_instructions=True,
        )

    # Retain upstream notices in the standalone artifact, even after minification.
    licenses = Path(HTML_ASSET_DIR) / 'licenses'
    notices = '\n\n'.join(
        f'{path.stem}\n{path.read_text(encoding="utf-8")}'
        for path in sorted(licenses.glob('*.txt'))
    )
    html_content += '\n<!-- Bundled library licenses\n' + notices.replace('--', '- -') + '\n-->\n'

    # Write the HTML file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
def main(argv=None):
    """Main function to process files and generate visualization."""
    parser = argparse.ArgumentParser(description='Generate an HTML-based genome coverage visualization')
    parser.add_argument('-r', '--results', '--result-dir', help='GOTTCHA2 result directory; discover matching files automatically')
    parser.add_argument('-p', '--prefix', help='Sample prefix to select when the result directory contains multiple samples')
    source = parser.add_mutually_exclusive_group()
    source.add_argument('-b', '--bam', help='Coordinate-sorted GOTTCHA2 BAM; generate coverage automatically')
    source.add_argument('-c', '--coverage', help='Existing samtools coverage TSV (skip coverage generation)')
    parser.add_argument('--reference', help='Signature reference FASTA (plain, gzip, or BGZF); enables haploid variant calling')
    parser.add_argument('-f', '--full',
                       help='Path to the full taxonomy profiling TSV with required columns')
    parser.add_argument('-o', '--output', help='Output HTML (default: <prefix>.coverage.html beside the results/BAM)')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Threads for variant calling (default: 1)')
    parser.add_argument('-e', '--external', action='store_true',
                       help='Deprecated compatibility flag; resources are always embedded')
    variants = parser.add_mutually_exclusive_group()
    variants.add_argument('--vcf', action='append', nargs='+', default=[],
                       help='Path(s) to bgzip/gzip-compressed .vcf.gz files. May be used more than once. Each VCF must have a matching .tbi index.')
    variants.add_argument('--no-variants', action='store_true', help='Generate coverage only, skipping VCF discovery and variant calling')
    parser.add_argument('--min-depth', type=int, default=5,
                       help='Minimum INFO/DP depth required to keep a VCF variant (default: 5)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error('--threads must be at least 1')
    if args.min_depth < 0:
        parser.error('--min-depth must be non-negative')

    # The renderer remains executable as a standalone script as well as a subcommand.
    if __package__:
        from .coverage_browser_workflow import resolve_inputs, prepare_inputs
        from .process_bam import write_taxid_network
    else:
        from coverage_browser_workflow import resolve_inputs, prepare_inputs
        from process_bam import write_taxid_network
    import pysam
    import logging

    if args.verbose:
        logging_level = logging.INFO
    else:
        logging_level = logging.WARNING

    logging.basicConfig(
        level=logging_level,
        format='[%(asctime)s] [%(levelname)s] [%(module)s] %(message)s',
        datefmt='%Y%m%d %H:%M:%S',
    )

    try:
        inputs = resolve_inputs(
            results=args.results, prefix=args.prefix, bam=args.bam, coverage=args.coverage,
            full=args.full, reference=args.reference, vcfs=_flatten_vcf_args(args.vcf),
            output=args.output, no_variants=args.no_variants,
        )

        _validate_vcf_inputs([str(path) for path in inputs.vcfs])
        inputs.output.parent.mkdir(parents=True, exist_ok=True)
        if not args.no_variants and not inputs.vcfs and inputs.reference is None:
            logging.warning('No reference FASTA or VCF found; generating coverage-only HTML.')
        with tempfile.TemporaryDirectory(prefix='.gottcha2-coverage-', dir=inputs.output.parent) as tmp:
            node_file = Path(args.results) / 'node.tsv'
            edge_file = Path(args.results) / 'edge.tsv' 
            if inputs.bam is not None:
                logging.info("Writing taxid network to node and edge files...")
                write_taxid_network(str(inputs.bam), str(node_file), str(edge_file))
            else:
                node_file = None
                edge_file = None
                logging.info("No BAM file provided; skipping taxid network generation...")
            logging.info("Preparing inputs for coverage analysis...")
            coverage, vcfs = prepare_inputs(inputs, tmp, threads=args.threads, no_variants=args.no_variants)
            logging.info("Parsing coverage file...")
            coverage_data, uniq_taxid_list = parse_coverage_file(coverage)
            logging.info("Parsing full genome file...")
            genome_data, species_data = parse_full_file(inputs.full, uniq_taxid_list)
            logging.info("Parsing VCF files...")
            variant_data, vcf_file_data = parse_vcf_files(
                [str(path) for path in vcfs], coverage_data, min_depth=args.min_depth
            )
            temporary_html = Path(tmp) / 'coverage.html'
            logging.info("Generating coverage HTML...")
            generate_html(
                coverage_data, genome_data, species_data, temporary_html,
                node_file=node_file, edge_file=edge_file, full_file=inputs.full,
                variant_data=variant_data, vcf_file_data=vcf_file_data,
                min_variant_depth=args.min_depth,
            )
            temporary_html.replace(inputs.output)
    except (OSError, ValueError, RuntimeError, pysam.SamtoolsError) as exc:
        parser.exit(1, f'Error: {exc}\n')

    logging.info("Coverage HTML generated: %s", inputs.output)

if __name__ == "__main__":
    main()
