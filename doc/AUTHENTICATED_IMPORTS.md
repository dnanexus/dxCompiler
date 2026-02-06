# Authenticated HTTP Imports for WDL

This document describes how to import WDL files from HTTP sources that require authentication.

## Overview

dxCompiler can import WDL files from URLs that require authentication, such as:
- Private GitHub repositories
- Private GitLab repositories (when configured)
- Internal corporate servers (when configured)

Authentication is provided via Bearer tokens set in environment variables.

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `WDL_IMPORT_TOKEN` | No | Bearer token for HTTP authentication |
| `WDL_IMPORT_TOKEN_DOMAINS` | No | Comma-separated list of domains to send token to |

## Default Allowed Domains

When `WDL_IMPORT_TOKEN_DOMAINS` is not set, tokens are only sent to:

- `github.com`
- `raw.githubusercontent.com`

This prevents accidentally leaking tokens to untrusted servers.

## Configuration Examples

### Basic GitHub Private Repository Access

```bash
# Generate a token at https://github.com/settings/tokens
# Required scope: repo (for private repositories)
export WDL_IMPORT_TOKEN="ghp_xxxxxxxxxxxxxxxxxxxx"

java -jar dxCompiler.jar compile workflow.wdl -project project-xxxx -folder /my/workflows/
```

### Multiple Private Sources

```bash
export WDL_IMPORT_TOKEN="your-token-here"
export WDL_IMPORT_TOKEN_DOMAINS="github.com,raw.githubusercontent.com,gitlab.com,internal.company.com"

java -jar dxCompiler.jar compile workflow.wdl -project project-xxxx -folder /my/workflows/
```

## WDL Import Syntax

### GitHub Raw Content URL

```wdl
import "https://raw.githubusercontent.com/owner/repo/branch/path/to/file.wdl"
```

### GitHub Blob URL (Not Recommended)

GitHub blob URLs (`github.com/owner/repo/blob/...`) do not return raw content.
Use raw.githubusercontent.com URLs instead.

## Error Messages

### 401 Unauthorized

```
HTTP 401 Unauthorized when accessing https://raw.githubusercontent.com/...
If this is a private repository, ensure WDL_IMPORT_TOKEN is set with a valid access token.
For GitHub: generate a token at https://github.com/settings/tokens with 'repo' scope.
```

**Solution:** Set the `WDL_IMPORT_TOKEN` environment variable with a valid token.

### 403 Forbidden

```
HTTP 403 Forbidden when accessing https://raw.githubusercontent.com/...
The token may be invalid or lack the required permissions.
```

**Solution:** The token may be invalid or lack the required permissions. For GitHub, ensure the token has the `repo` scope.

## Security Considerations

1. **Token Scope**: Only grant the minimum required permissions to your token
2. **Domain Allowlist**: Tokens are only sent to explicitly allowed domains
3. **No Logging**: Token values are never logged; only usage is traced
4. **HTTPS Recommended**: Always use HTTPS URLs for private imports

## Debugging

Enable trace logging to see when authenticated imports are used:

```bash
java -jar dxCompiler.jar compile workflow.wdl -verbose -verboseKey FileSourceResolver
```

Look for log entries like:
```
[TRACE] WDL_IMPORT_TOKEN found; authenticated HTTP imports enabled for domains: github.com, raw.githubusercontent.com
[TRACE] Using authenticated HTTP for import from: raw.githubusercontent.com
```

## Troubleshooting

### Token not being sent

1. Verify `WDL_IMPORT_TOKEN` is set: `echo $WDL_IMPORT_TOKEN`
2. Check if the domain is in the allowed list
3. Enable verbose logging to see authentication attempts

### Token rejected

1. Verify the token is still valid (not expired)
2. Check token permissions/scopes
3. For GitHub, ensure the token has access to the specific repository

### Public imports stopped working

The authenticated HTTP protocol is backward compatible. If no token is set, it works like the standard HTTP protocol. Verify no token is set if you're testing public access.

## How It Works

dxCompiler uses a custom `AuthenticatedHttpFileAccessProtocol` that:

1. Checks if `WDL_IMPORT_TOKEN` is set
2. For each HTTP import, checks if the domain is in the allowed list
3. If both conditions are met, adds an `Authorization: Bearer <token>` header to the request
4. If either condition is not met, the request proceeds without authentication (backward compatible)

This ensures that:
- Tokens are never sent to untrusted domains
- Existing workflows continue to work without modification
- Authentication failures produce clear, actionable error messages
